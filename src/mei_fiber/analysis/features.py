"""Per-fiber feature-fraction export for ECDF-style summaries."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Any

import numpy as np

from mei_fiber.db import FiberDatabase

_ALL_LAYERS = ["nucleosomes", "5mC", "5hmC", "6mA", "msp"]


def _slice(chrom_data: dict[str, Any], layer: str, fiber_int_id: int) -> tuple[int, int] | None:
    return chrom_data.get(f"{layer}_idx", {}).get(int(fiber_int_id))


def _interval_coverage(
    starts: np.ndarray,
    ends: np.ndarray,
    region_start: int,
    region_end: int,
) -> int:
    left = np.maximum(starts.astype(np.int64), region_start)
    right = np.minimum(ends.astype(np.int64), region_end)
    return int(np.maximum(right - left, 0).sum())


def _safe_fraction(numerator: float, denominator: float) -> float:
    return float(numerator / denominator) if denominator > 0 else 0.0


def _fiber_feature_row(
    chrom_data: dict[str, Any],
    sample: str,
    chrom: str,
    fiber_int_id: int,
    fiber_start: int,
    fiber_end: int,
) -> dict[str, Any]:
    fiber_len = max(int(fiber_end) - int(fiber_start), 0)
    row: dict[str, Any] = {
        "sample": sample,
        "chrom": chrom,
        "fiber_id": chrom_data["id_table"][int(fiber_int_id)].decode(),
        "fiber_start": int(fiber_start),
        "fiber_end": int(fiber_end),
        "fiber_length": fiber_len,
    }

    nuc_bp = 0
    nuc_slice = _slice(chrom_data, "nucleosomes", fiber_int_id)
    if nuc_slice and "nucleosomes" in chrom_data:
        s, e = nuc_slice
        nuc = chrom_data["nucleosomes"]
        nuc_bp = _interval_coverage(nuc["starts"][s:e], nuc["ends"][s:e], fiber_start, fiber_end)
    row["frac_nuc_coverage"] = _safe_fraction(nuc_bp, fiber_len)

    msp_bp = 0
    msp_slice = _slice(chrom_data, "msp", fiber_int_id)
    if msp_slice and "msp" in chrom_data:
        s, e = msp_slice
        msp = chrom_data["msp"]
        msp_bp = _interval_coverage(msp["starts"][s:e], msp["ends"][s:e], fiber_start, fiber_end)
    row["frac_msp_coverage"] = _safe_fraction(msp_bp, fiber_len)

    for layer, count_col, fraction_col in [
        ("5mC", "n_5mC", "frac_5mC_coverage"),
        ("5hmC", "n_5hmC", "frac_5hmC_coverage"),
        ("6mA", "n_6mA", "frac_6mA_coverage"),
    ]:
        sl = _slice(chrom_data, layer, fiber_int_id)
        n_calls = 0
        if sl and layer in chrom_data:
            s, e = sl
            positions = chrom_data[layer]["positions"][s:e]
            n_calls = int(((positions >= fiber_start) & (positions < fiber_end)).sum())
        row[count_col] = n_calls
        row[fraction_col] = _safe_fraction(n_calls, fiber_len)

    for layer, output_col in [
        ("5mC", "frac_5mC_among_cpg"),
        ("5hmC", "frac_5hmC_among_cpg"),
    ]:
        sl = _slice(chrom_data, layer, fiber_int_id)
        value = np.nan
        if sl and layer in chrom_data and "is_methylated" in chrom_data[layer]:
            s, e = sl
            calls = chrom_data[layer]["is_methylated"][s:e]
            positions = chrom_data[layer]["positions"][s:e]
            keep = (positions >= fiber_start) & (positions < fiber_end)
            if keep.any():
                value = float(calls[keep].mean())
        row[output_col] = value

    # The database stores called 6mA positions, not all adenines. This column is kept
    # for compatibility with the older ECDF notebook, but it is a density-style proxy.
    row["frac_6mA_among_A"] = row["frac_6mA_coverage"]
    return row


def export_global_feature_fractions(
    db_path: Path,
    out_csv: Path,
    samples: list[str] | None = None,
    max_fibers_per_chrom: int | None = None,
) -> Path:
    """Export global per-fiber feature fractions for ECDF plotting.

    Fractions are computed across each fiber's full span. Nucleosome/MSP fractions are
    covered base pairs divided by fiber length. Modification coverage fractions are
    call counts divided by fiber length; 5mC/5hmC methylation fractions use the stored
    binary calls.
    """
    out_csv = Path(out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)

    fieldnames = [
        "sample",
        "chrom",
        "fiber_id",
        "fiber_start",
        "fiber_end",
        "fiber_length",
        "frac_5mC_coverage",
        "frac_5hmC_coverage",
        "frac_6mA_coverage",
        "frac_msp_coverage",
        "frac_nuc_coverage",
        "frac_5mC_among_cpg",
        "frac_5hmC_among_cpg",
        "frac_6mA_among_A",
        "n_5mC",
        "n_5hmC",
        "n_6mA",
    ]
    with FiberDatabase(db_path) as db, out_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        selected_samples = samples or db.samples
        for sample in selected_samples:
            if sample not in db.samples:
                raise ValueError(f"Unknown sample {sample!r}; available: {db.samples}")
            for chrom in db.get_chromosomes(sample):
                chrom_data = db.load_chromosome_data(sample, chrom, layers=_ALL_LAYERS)
                if chrom_data is None:
                    continue
                limit = len(chrom_data["fm_ids"])
                if max_fibers_per_chrom is not None:
                    limit = min(limit, max_fibers_per_chrom)
                for fid, start, end in zip(
                    chrom_data["fm_ids"][:limit],
                    chrom_data["fm_starts"][:limit],
                    chrom_data["fm_ends"][:limit],
                    strict=False,
                ):
                    writer.writerow(
                        _fiber_feature_row(chrom_data, sample, chrom, int(fid), int(start), int(end))
                    )
    return out_csv


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--samples", nargs="+")
    parser.add_argument("--max-fibers-per-chrom", type=int)
    args = parser.parse_args()

    out = export_global_feature_fractions(
        args.db,
        args.out,
        samples=args.samples,
        max_fibers_per_chrom=args.max_fibers_per_chrom,
    )
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
