"""PacBio FIRE co-accessibility helpers.

This module ports the first, testable pieces of the legacy FIRE co-accessibility
workflow into PACKAGE:

1. Build a ``Cov.bed``-style table from a PACKAGE HDF5 database that contains the
   ``fire_accessibility`` layer.
2. Convert that ``Cov.bed`` table into the enhancer-by-fiber JSON object used by the
   original ranking script.
3. Rank constituent FIRE element pairs using the legacy co-accessibility score.

The generated ``Cov.bed`` has the same nine columns as the old
``bedtools intersect + awk`` pipeline:

``element_chr element_start element_end stitch_chr stitch_start stitch_end fiber_id fire_score overlap_bp``
"""

from __future__ import annotations

import argparse
import json
import math
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


@dataclass(frozen=True)
class CoaccessPairScore:
    """Rankable co-accessibility score for one pair of constituent FIRE elements."""

    score: float
    stitched_id: str
    element_a: str
    element_b: str
    rank_counts: tuple[int, int, int, int]
    distance_bp: int


def _powlaw(x: np.ndarray | float, a: float, b: float, c: float) -> np.ndarray | float:
    return a / np.power(x, b) + c


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
    text = f"{float(score):.9g}"
    if "." in text and "e" not in text.lower():
        text = text.rstrip("0").rstrip(".")
    return text or "0"


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


def _element_distance_bp(element_a: str, element_b: str) -> int:
    chrom_a, coords_a = element_a.split(":")
    chrom_b, coords_b = element_b.split(":")
    if chrom_a != chrom_b:
        return 0
    start_a, end_a = [int(x) for x in coords_a.split("-")]
    start_b, end_b = [int(x) for x in coords_b.split("-")]
    if end_a > end_b:
        return start_a - end_b
    return start_b - end_a


def _raw_pair_score(
    fibers_a: list[float | None],
    fibers_b: list[float | None],
    *,
    threshold: float,
) -> tuple[float | None, tuple[int, int, int, int]]:
    """Compute the legacy modified odds-ratio-like score for one element pair.

    Counts match ``fireAB_ranking.py``:

    - ``a``: both elements accessible on the same fiber (score < threshold)
    - ``b``/``c``: discordant accessibility
    - ``d``: both present but not accessible

    The old script keeps pairs with at least 10 shared fibers, at least one
    discordant fiber, and at least one concordant accessible fiber.
    """
    a = b = c = d = 0
    for score_a, score_b in zip(fibers_a, fibers_b, strict=False):
        if score_a is None or score_b is None:
            continue
        accessible_a = score_a < threshold
        accessible_b = score_b < threshold
        if accessible_a and accessible_b:
            a += 1
        elif accessible_a and not accessible_b:
            c += 1
        elif not accessible_a and accessible_b:
            b += 1
        else:
            d += 1

    if a + b + c + d < 10 or b + c == 0 or a == 0:
        return None, (a, b, c, d)

    discordant_mean = (b + c) / 2
    return (a**2) / (discordant_mean**2), (a, b, c, d)


def _iter_pair_scores(
    matrix: list[dict[str, Any]],
    *,
    threshold: float,
    distance_params: tuple[float, float, float] | None = None,
) -> list[CoaccessPairScore]:
    scores: list[CoaccessPairScore] = []
    for stitched in matrix:
        stitched_id = stitched["seId"]
        enhancers = stitched.get("enhs", [])
        if len(enhancers) < 2:
            continue
        for i, enh_a in enumerate(enhancers):
            for j in range(i + 1, len(enhancers)):
                enh_b = enhancers[j]
                raw_score, counts = _raw_pair_score(
                    enh_a["fibers"],
                    enh_b["fibers"],
                    threshold=threshold,
                )
                if raw_score is None:
                    continue
                distance = _element_distance_bp(enh_a["enhId"], enh_b["enhId"])
                score = raw_score
                if distance_params is not None:
                    if distance > 50:
                        expected = float(_powlaw(distance, *distance_params))
                        score = raw_score / expected if expected > 0 else 0
                    else:
                        score = 0
                scores.append(
                    CoaccessPairScore(
                        score=float(score),
                        stitched_id=stitched_id,
                        element_a=enh_a["enhId"],
                        element_b=enh_b["enhId"],
                        rank_counts=counts,
                        distance_bp=distance,
                    )
                )
    return scores


def _calculate_elbow_cutoff(values: list[float]) -> float:
    """Return the same y-axis cutoff used by the legacy ranking plots."""
    if not values:
        return math.inf
    arr = np.sort(np.asarray(values, dtype=float))
    arr[arr < 0] = 0
    if len(arr) == 1 or np.isclose(arr.max(), arr.min()):
        return float(arr[-1])

    slope = (float(arr.max()) - float(arr.min())) / len(arr)
    x_pts = np.arange(1, len(arr) + 1)
    below_counts = []
    for idx, y_pt in enumerate(arr):
        b = y_pt - (slope * idx)
        below_counts.append(int(np.sum(arr <= (x_pts * slope + b))))
    x_pt = int(np.floor(np.argmin(below_counts)) + 1)
    return float(arr[x_pt - 1])


def _write_elbow_plot(values: list[float], out_svg: Path, title: str) -> None:
    if not values:
        return
    try:
        from matplotlib import pyplot as plt
    except ImportError:
        return

    arr = np.sort(np.asarray(values, dtype=float))
    cutoff = _calculate_elbow_cutoff(values)
    x_pt = int(np.searchsorted(arr, cutoff, side="left")) + 1
    out_svg.parent.mkdir(parents=True, exist_ok=True)

    plt.figure(figsize=(10, 6))
    plt.plot(np.arange(1, len(arr) + 1), arr, "b-")
    plt.scatter(x_pt, cutoff, marker="o", color="red")
    median = float(np.median(arr)) if len(arr) else float("nan")
    mean = float(np.mean(arr)) if len(arr) else float("nan")
    median_fold = cutoff / median if median else float("inf")
    mean_fold = cutoff / mean if mean else float("inf")
    plt.title(
        f"{title}\nx={x_pt}\ny={cutoff:.3f}\n"
        f"Fold over Median={median_fold:.3f}x\nFold over Mean={mean_fold:.3f}x"
    )
    plt.xlabel("Elements ranked by modified odds ratio")
    plt.ylabel("Modified odds ratio")
    plt.tight_layout()
    plt.savefig(out_svg)
    plt.close()


def _fit_distance_correction(
    pair_scores: list[CoaccessPairScore],
    *,
    window_size: int = 500,
) -> tuple[float, float, float] | None:
    """Fit the legacy power-law distance correction from uncorrected pair scores."""
    rows = sorted(
        [(s.distance_bp, s.score) for s in pair_scores if s.distance_bp > 0 and s.score > 0],
        key=lambda x: x[0],
    )
    if len(rows) < window_size:
        return None

    x: list[float] = []
    y: list[float] = []
    for i in range(0, len(rows), window_size):
        chunk = rows[i : i + window_size]
        if len(chunk) < window_size:
            continue
        distances, scores = zip(*chunk, strict=False)
        x.append(float(np.median(distances)))
        y.append(float(np.median(scores)))

    if len(x) < 3:
        return None

    x_arr = np.asarray(x, dtype=float)
    y_arr = np.asarray(y, dtype=float)
    try:
        from scipy.optimize import curve_fit

        params, _ = curve_fit(_powlaw, x_arr, y_arr, maxfev=10_000_000)
        return tuple(float(v) for v in params)
    except Exception:
        # Fallback for environments without scipy. This is not byte-for-byte legacy,
        # but preserves the same monotonic distance-correction idea.
        baseline = max(float(np.percentile(y_arr, 10)), 1e-9)
        shifted = np.maximum(y_arr - baseline, 1e-9)
        slope, intercept = np.polyfit(np.log(x_arr), np.log(shifted), 1)
        return (float(np.exp(intercept)), float(-slope), baseline)


def rank_coaccessibility_object(
    object_json: Path,
    outdir: Path,
    *,
    threshold: float = 0.10,
    distance_correct: bool = True,
    write_plots: bool = True,
) -> dict[str, Any]:
    """Rank constituent FIRE element pairs from a co-accessibility object.

    This ports the ranking endpoint of the legacy ``fireAB_ranking.py`` workflow.
    It writes:

    - ``ce_rank.txt``: ranked constituent-element pairs.
    - ``cluster_rank.txt``: one summary row per stitched region, using the maximum
      pair score in that region. The legacy script had clustering disabled, so this
      corresponds to its default behavior.
    - ``ce_pairs_ranked.svg`` and ``clusters_ranked.svg`` when matplotlib is
      available and ``write_plots`` is true.
    """
    with Path(object_json).open() as handle:
        matrix = json.load(handle)

    raw_scores = _iter_pair_scores(matrix, threshold=threshold)
    distance_params = _fit_distance_correction(raw_scores) if distance_correct else None
    ranked_pairs = _iter_pair_scores(
        matrix,
        threshold=threshold,
        distance_params=distance_params,
    )
    ranked_pairs.sort(key=lambda s: s.score, reverse=True)

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    pair_cutoff = _calculate_elbow_cutoff([s.score for s in ranked_pairs])
    pair_lines = []
    for rank, score in enumerate(ranked_pairs, start=1):
        is_super = "Super" if score.score >= pair_cutoff else ""
        pair_lines.append(
            "\t".join(
                [
                    score.element_a,
                    score.element_b,
                    score.stitched_id,
                    str(score.score),
                    str(rank),
                    is_super,
                ]
            )
            + "\n"
        )
    (outdir / "ce_rank.txt").write_text("".join(pair_lines))

    by_stitched: dict[str, list[CoaccessPairScore]] = defaultdict(list)
    for score in ranked_pairs:
        by_stitched[score.stitched_id].append(score)

    enh_by_stitched = {
        row["seId"]: [enh["enhId"] for enh in row.get("enhs", [])]
        for row in matrix
    }
    cluster_rows: list[tuple[float, str, list[str], tuple[int, int, int, int]]] = []
    for stitched_id, scores in by_stitched.items():
        best = max(scores, key=lambda s: s.score)
        cluster_rows.append(
            (
                best.score,
                stitched_id,
                enh_by_stitched.get(stitched_id, []),
                best.rank_counts,
            )
        )
    cluster_rows.sort(key=lambda row: row[0], reverse=True)

    cluster_cutoff = _calculate_elbow_cutoff([row[0] for row in cluster_rows])
    cluster_lines = []
    for rank, (score, stitched_id, enh_ids, counts) in enumerate(cluster_rows, start=1):
        is_super = "Super" if score >= cluster_cutoff else ""
        ce_string = ", ".join(enh_ids)
        cluster_lines.append(
            "\t".join(
                [
                    ce_string,
                    stitched_id,
                    str(len(enh_ids)),
                    str(score),
                    str(rank),
                    is_super,
                    *(str(v) for v in counts),
                ]
            )
            + "\n"
        )
    (outdir / "cluster_rank.txt").write_text("".join(cluster_lines))

    if write_plots:
        _write_elbow_plot(
            [s.score for s in ranked_pairs],
            outdir / "ce_pairs_ranked.svg",
            "Constituent FIRE pairs",
        )
        _write_elbow_plot(
            [row[0] for row in cluster_rows],
            outdir / "clusters_ranked.svg",
            "Stitched FIRE regions",
        )

    return {
        "object_json": str(object_json),
        "outdir": str(outdir),
        "threshold": threshold,
        "distance_correct": distance_correct,
        "distance_params": distance_params,
        "n_ranked_pairs": len(ranked_pairs),
        "n_ranked_clusters": len(cluster_rows),
        "pair_cutoff": pair_cutoff,
        "cluster_cutoff": cluster_cutoff,
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

    rank = subparsers.add_parser("rank", help="Rank co-accessible FIRE element pairs")
    rank.add_argument("--object", type=Path, required=True)
    rank.add_argument("--outdir", type=Path, required=True)
    rank.add_argument("--threshold", type=float, default=0.10)
    rank.add_argument("--no-distance-correct", action="store_true")
    rank.add_argument("--no-plots", action="store_true")

    args = parser.parse_args()
    if args.command == "cov":
        summary = export_coaccessibility_cov(
            args.db,
            args.stitched,
            args.peaks,
            args.out,
            sample=args.sample,
        )
    elif args.command == "object":
        summary = cov_to_object(args.cov, args.out)
    else:
        summary = rank_coaccessibility_object(
            args.object,
            args.outdir,
            threshold=args.threshold,
            distance_correct=not args.no_distance_correct,
            write_plots=not args.no_plots,
        )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
