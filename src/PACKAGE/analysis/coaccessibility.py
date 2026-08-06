"""PacBio FIRE co-accessibility helpers.

This module ports the first, testable pieces of the legacy FIRE co-accessibility
workflow into PACKAGE:

1. Build a ``Cov.bed``-style table from a PACKAGE HDF5 database that contains the
   ``fire_accessibility`` layer.
2. Convert that ``Cov.bed`` table into the enhancer-by-fiber JSON object used by the
   original ranking script.

The generated ``Cov.bed`` has the same nine columns as the old
``bedtools intersect + awk`` pipeline:

``element_chr element_start element_end stitch_chr stitch_start stitch_end fiber_id fire_score overlap_bp``
"""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from PACKAGE.db import FiberDatabase


@dataclass(frozen=True, order=True)
class BedInterval:
    """Simple 0-based half-open BED interval."""

    chrom: str
    start: int
    end: int

    @property
    def region_id(self) -> str:
        return f"{self.chrom}:{self.start}-{self.end}"


@dataclass(frozen=True)
class CoaccessRegion:
    """Constituent FIRE element nested inside a stitched FIRE region."""

    element: BedInterval
    stitched: BedInterval


def _read_bed3(path: Path) -> list[BedInterval]:
    intervals: list[BedInterval] = []
    with Path(path).open() as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 3:
                raise ValueError(f"Expected BED3+ row at {path}:{line_number}")
            try:
                start = int(cols[1])
                end = int(cols[2])
            except ValueError as exc:
                raise ValueError(f"Invalid BED coordinates at {path}:{line_number}") from exc
            if end < start:
                raise ValueError(f"Invalid BED interval with end < start at {path}:{line_number}")
            intervals.append(BedInterval(cols[0], start, end))
    return sorted(intervals)


def _group_by_chrom(intervals: list[BedInterval]) -> dict[str, list[BedInterval]]:
    grouped: dict[str, list[BedInterval]] = defaultdict(list)
    for interval in intervals:
        grouped[interval.chrom].append(interval)
    return {chrom: sorted(rows) for chrom, rows in grouped.items()}


def find_contained_peak_pairs(
    stitched_bed: Path,
    peaks_bed: Path,
) -> list[CoaccessRegion]:
    """Return FIRE peak/stitch pairs matching the old first bedtools intersect.

    The legacy command was:

    ``bedtools intersect -sorted -F 1 -a FIRE_stitched.bed -b FIRE_peaks_intergenic.bed``

    With stitched intervals as ``-a`` and peaks as ``-b``, ``-F 1`` means the full
    peak must be contained in the stitched interval. The legacy awk then printed
    the peak coordinates first and stitched coordinates second.
    """
    stitched_by_chrom = _group_by_chrom(_read_bed3(stitched_bed))
    peaks_by_chrom = _group_by_chrom(_read_bed3(peaks_bed))

    pairs: list[CoaccessRegion] = []
    for chrom, peaks in peaks_by_chrom.items():
        stitched = stitched_by_chrom.get(chrom, [])
        active: list[BedInterval] = []
        stitch_i = 0
        for peak in peaks:
            while stitch_i < len(stitched) and stitched[stitch_i].start <= peak.start:
                active.append(stitched[stitch_i])
                stitch_i += 1
            active = [region for region in active if region.end > peak.start]
            for region in active:
                if region.start <= peak.start and region.end >= peak.end:
                    pairs.append(CoaccessRegion(element=peak, stitched=region))
    return sorted(pairs, key=lambda p: (p.element.chrom, p.element.start, p.element.end, p.stitched.start, p.stitched.end))


def _format_score(score: float) -> str:
    """Format FIRE scores close to the legacy text output."""
    return f"{float(score):.9g}"


def _write_cov_for_chrom(
    db: FiberDatabase,
    sample: str,
    chrom: str,
    pairs: list[CoaccessRegion],
    handle: Any,
) -> int:
    chrom_data = db.load_chromosome_data(sample, chrom, layers=["fire_accessibility"])
    if chrom_data is None or "fire_accessibility" not in chrom_data:
        return 0

    fire = chrom_data["fire_accessibility"]
    if len(fire["starts"]) == 0:
        return 0

    fiber_ids = fire["fiber_int_ids"]
    starts = fire["starts"]
    ends = fire["ends"]
    scores = fire["scores"]
    id_table = chrom_data["id_table"]

    order = np.lexsort((fiber_ids, ends, starts))
    starts = starts[order]
    ends = ends[order]
    scores = scores[order]
    fiber_ids = fiber_ids[order]

    cursor = 0
    active: list[int] = []
    rows = 0
    for pair in pairs:
        peak = pair.element
        stitched = pair.stitched
        while cursor < len(starts) and int(starts[cursor]) < peak.end:
            active.append(cursor)
            cursor += 1
        active = [i for i in active if int(ends[i]) > peak.start]
        for i in active:
            fire_start = int(starts[i])
            fire_end = int(ends[i])
            if fire_start >= peak.end or fire_end <= peak.start:
                continue
            overlap = min(fire_end, peak.end) - max(fire_start, peak.start)
            if overlap <= 0:
                continue
            fiber_id = id_table[int(fiber_ids[i])].decode()
            handle.write(
                "\t".join(
                    [
                        peak.chrom,
                        str(peak.start),
                        str(peak.end),
                        stitched.chrom,
                        str(stitched.start),
                        str(stitched.end),
                        fiber_id,
                        _format_score(float(scores[i])),
                        str(overlap),
                    ]
                )
                + "\n"
            )
            rows += 1
    return rows


def export_coaccessibility_cov(
    db_path: Path,
    stitched_bed: Path,
    peaks_bed: Path,
    out_bed: Path,
    *,
    sample: str | None = None,
) -> dict[str, int | str]:
    """Export legacy-compatible ``Cov.bed`` from PACKAGE HDF5 fire calls.

    Args:
        db_path: PACKAGE HDF5 database containing ``fire_accessibility``.
        stitched_bed: ``FIRE_stitched.bed`` from the legacy stitcher.
        peaks_bed: ``FIRE_peaks_intergenic.bed`` or equivalent constituent FIRE peaks.
        out_bed: Output path for the nine-column ``Cov.bed`` table.
        sample: Sample name. Defaults to the first sample in the database.

    Returns:
        Summary dict with region-pair and output-row counts.
    """
    pairs = find_contained_peak_pairs(stitched_bed, peaks_bed)
    pairs_by_chrom: dict[str, list[CoaccessRegion]] = defaultdict(list)
    for pair in pairs:
        pairs_by_chrom[pair.element.chrom].append(pair)

    out_bed = Path(out_bed)
    out_bed.parent.mkdir(parents=True, exist_ok=True)
    rows = 0
    with FiberDatabase(db_path) as db, out_bed.open("w") as handle:
        selected_sample = sample or db.samples[0]
        if selected_sample not in db.samples:
            raise ValueError(f"Unknown sample {selected_sample!r}; available: {db.samples}")
        for chrom in sorted(pairs_by_chrom):
            rows += _write_cov_for_chrom(db, selected_sample, chrom, pairs_by_chrom[chrom], handle)

    return {
        "db": str(db_path),
        "sample": sample or "",
        "stitched_bed": str(stitched_bed),
        "peaks_bed": str(peaks_bed),
        "out_bed": str(out_bed),
        "contained_peak_pairs": len(pairs),
        "cov_rows": rows,
    }


def cov_to_object(cov_bed: Path, out_json: Path) -> dict[str, int | str]:
    """Convert ``Cov.bed`` to the legacy enhancer-by-fiber JSON object.

    This mirrors ``objPrep.py`` but compares FIRE scores numerically when selecting
    the best overlapping score for a fiber/element pair.
    """
    se_fire_score: dict[str, dict[str, Any]] = {}
    with Path(cov_bed).open() as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            cells = line.rstrip("\n").split("\t")
            if len(cells) != 9:
                raise ValueError(f"Expected 9 columns in Cov.bed at {cov_bed}:{line_number}")
            enh_id = f"{cells[0]}:{cells[1]}-{cells[2]}"
            super_id = f"{cells[3]}:{cells[4]}-{cells[5]}"
            fiber_id = cells[6]
            fire_score = float(cells[7])
            feature_length = int(cells[8])

            if super_id not in se_fire_score:
                se_fire_score[super_id] = {"fiberList": []}
            if enh_id not in se_fire_score[super_id]:
                se_fire_score[super_id][enh_id] = {}
            if fiber_id not in se_fire_score[super_id]["fiberList"]:
                se_fire_score[super_id]["fiberList"].append(fiber_id)

            fiber_entry = se_fire_score[super_id][enh_id].get(fiber_id)
            if fiber_entry is None:
                se_fire_score[super_id][enh_id][fiber_id] = {
                    "best_score": fire_score,
                    "length": feature_length,
                }
            else:
                if fire_score < fiber_entry["best_score"]:
                    fiber_entry["best_score"] = fire_score
                fiber_entry["length"] += feature_length

    matrix: list[dict[str, Any]] = []
    for se_id, se in se_fire_score.items():
        se_mat = {"seId": se_id, "enhs": []}
        for enh_id, enh in se.items():
            if enh_id == "fiberList":
                continue
            enh_start, enh_end = enh_id.split(":")[1].split("-")
            enh_size = int(enh_end) - int(enh_start)
            enh_column: dict[str, Any] = {"enhId": enh_id, "fibers": []}
            for fiber in se["fiberList"]:
                fiber_entry = enh.get(fiber)
                if fiber_entry is None:
                    enh_column["fibers"].append(None)
                    continue
                fire_score = float(fiber_entry["best_score"])
                if fire_score < 0.1 or fiber_entry["length"] == enh_size:
                    enh_column["fibers"].append(fire_score)
                else:
                    enh_column["fibers"].append(None)
            se_mat["enhs"].append(enh_column)
        matrix.append(se_mat)

    out_json = Path(out_json)
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(matrix))
    return {
        "cov_bed": str(cov_bed),
        "out_json": str(out_json),
        "stitched_regions": len(matrix),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    cov = subparsers.add_parser("cov", help="Export legacy-compatible Cov.bed")
    cov.add_argument("--db", type=Path, required=True)
    cov.add_argument("--sample")
    cov.add_argument("--stitched", type=Path, required=True)
    cov.add_argument("--peaks", type=Path, required=True)
    cov.add_argument("--out", type=Path, required=True)

    obj = subparsers.add_parser("object", help="Convert Cov.bed to legacy JSON object")
    obj.add_argument("--cov", type=Path, required=True)
    obj.add_argument("--out", type=Path, required=True)

    args = parser.parse_args()
    if args.command == "cov":
        summary = export_coaccessibility_cov(
            args.db,
            args.stitched,
            args.peaks,
            args.out,
            sample=args.sample,
        )
    else:
        summary = cov_to_object(args.cov, args.out)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
