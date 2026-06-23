"""Annotation-centered methylation and nucleosome matrix export."""

from __future__ import annotations

import argparse
import csv
import warnings
from pathlib import Path
from typing import Any

import numpy as np

from PACKAGE.db import FiberDatabase

BODY_BINS_CONFIG: dict[str, int] = {
    "Promoter": 2500,
    "Bivalent_Promoter": 2500,
    "PRC_Promoter": 2500,
    "CGI": 1000,
    "DMV": 7000,
    "Typical_Enhancer": 1000,
    "Super_Enhancer": 10000,
    "Gene_Body": 30000,
    "Bivalent_Genebody": 30000,
    "PRC_Genebody": 30000,
    "FIRE_d0_Peak": 500,
    "FIRE_d4_Peak": 500,
}
DEFAULT_BODY_BINS = 2500
METAPLOT_BODY_BINS = 200
FLANK_BINS = 100
FLANK_BP = 2000
MAX_REGION_BP = 50_000
MIN_FIBER_COVERAGE = 0.8
MIN_FIBERS = 10
MET_SMOOTHING_SIGMA = 50


def resize_vector(values: np.ndarray, target_bins: int) -> np.ndarray:
    """Resize one vector to a fixed number of bins by averaging source intervals."""
    n = len(values)
    if n == 0:
        return np.full(target_bins, np.nan)
    if n == target_bins:
        return values.astype(np.float64)

    result = np.full(target_bins, np.nan, dtype=np.float64)
    edges = np.linspace(0, n, target_bins + 1)
    values = values.astype(np.float64)
    for i in range(target_bins):
        lo = int(np.floor(edges[i]))
        hi = int(np.ceil(edges[i + 1]))
        hi = min(max(hi, lo + 1), n)
        chunk = values[lo:hi]
        if len(chunk) == 0:
            continue
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            result[i] = np.nanmean(chunk)
    return result


def _gaussian_kernel(sigma: float) -> np.ndarray:
    radius = max(1, int(4 * sigma))
    x = np.arange(-radius, radius + 1, dtype=np.float64)
    kernel = np.exp(-(x * x) / (2 * sigma * sigma))
    return kernel / kernel.sum()


def _same_length_convolve(values: np.ndarray, kernel: np.ndarray) -> np.ndarray:
    """Convolve while preserving the length of ``values`` even for short vectors."""
    smoothed = np.convolve(values, kernel, mode="same")
    target = len(values)
    if len(smoothed) == target:
        return smoothed
    if len(smoothed) > target:
        extra = len(smoothed) - target
        left = extra // 2
        return smoothed[left:left + target]
    pad_left = (target - len(smoothed)) // 2
    pad_right = target - len(smoothed) - pad_left
    return np.pad(smoothed, (pad_left, pad_right), mode="edge")


def smooth_methylation(values: np.ndarray, sigma: float = MET_SMOOTHING_SIGMA) -> np.ndarray:
    """Gaussian-smooth a methylation vector while preserving all-NaN positions."""
    values = values.astype(np.float64)
    if len(values) == 0 or np.isnan(values).all() or sigma <= 0:
        return values.copy()
    nan_mask = np.isnan(values)
    filled = values.copy()
    if nan_mask.any():
        valid = np.where(~nan_mask)[0]
        filled[nan_mask] = np.interp(np.where(nan_mask)[0], valid, filled[valid])
    smoothed = _same_length_convolve(filled, _gaussian_kernel(sigma))
    smoothed[nan_mask] = np.nan
    return smoothed


def _extract_nuc_occupancy(
    chrom_data: dict[str, Any],
    fiber_int_id: int,
    start: int,
    length: int,
) -> np.ndarray:
    occupancy = np.zeros(length, dtype=np.uint8)
    sl = chrom_data.get("nucleosomes_idx", {}).get(int(fiber_int_id))
    if sl is None or "nucleosomes" not in chrom_data:
        return occupancy
    s, e = sl
    data = chrom_data["nucleosomes"]
    keep = (data["ends"][s:e] > start) & (data["starts"][s:e] < start + length)
    for nuc_start, nuc_end in zip(data["starts"][s:e][keep], data["ends"][s:e][keep], strict=False):
        left = max(0, int(nuc_start) - start)
        right = min(length, int(nuc_end) - start)
        if right > left:
            occupancy[left:right] = 1
    return occupancy


def _extract_5mc(
    chrom_data: dict[str, Any],
    fiber_int_id: int,
    start: int,
    length: int,
    methylation_threshold: float,
) -> np.ndarray:
    methylation = np.full(length, np.nan, dtype=np.float64)
    sl = chrom_data.get("5mC_idx", {}).get(int(fiber_int_id))
    if sl is None or "5mC" not in chrom_data:
        return methylation
    s, e = sl
    data = chrom_data["5mC"]
    positions = data["positions"][s:e]
    keep = (positions >= start) & (positions < start + length)
    if not keep.any():
        return methylation
    rel_pos = positions[keep] - start
    if "is_methylated" in data:
        calls = data["is_methylated"][s:e][keep].astype(np.float64)
    elif "probabilities" in data:
        calls = (data["probabilities"][s:e][keep] >= methylation_threshold).astype(np.float64)
    else:
        return methylation
    methylation[rel_pos] = calls
    return methylation


def _nanmean_matrix(rows: list[np.ndarray]) -> np.ndarray:
    matrix = np.array(rows, dtype=np.float64)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        return np.nanmean(matrix, axis=0)


def _write_region_summary(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def export_annotation_matrices(
    db_path: Path,
    outdir: Path,
    annotations: list[str],
    samples: list[str] | None = None,
    max_regions: int | None = None,
    min_fiber_coverage: float = MIN_FIBER_COVERAGE,
    min_fibers: int = MIN_FIBERS,
    max_region_bp: int = MAX_REGION_BP,
    flank_bp: int = FLANK_BP,
) -> list[dict[str, Any]]:
    """Export normalized 5mC/nucleosome matrices for annotation-centered plots.

    For each annotation/sample pair, writes:
    ``<annotation>_<sample>_nuc_body.npz``, ``_met_body.npz``,
    ``_nuc_metaplot.npz``, and ``_met_metaplot.npz``.
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    all_summary: list[dict[str, Any]] = []

    with FiberDatabase(db_path) as db:
        selected_samples = samples or db.samples
        for sample in selected_samples:
            if sample not in db.samples:
                raise ValueError(f"Unknown sample {sample!r}; available: {db.samples}")

        for annotation in annotations:
            regions = db.get_annotation_regions_with_ids(annotation)
            if not regions:
                raise ValueError(f"Unknown or empty annotation {annotation!r}")
            if max_regions is not None:
                regions = regions[:max_regions]
            grouped = db.group_regions_by_chrom(regions)
            body_bins = BODY_BINS_CONFIG.get(annotation, DEFAULT_BODY_BINS)

            for sample in selected_samples:
                nuc_body: list[np.ndarray] = []
                met_body: list[np.ndarray] = []
                nuc_meta: list[np.ndarray] = []
                met_meta: list[np.ndarray] = []

                for chrom, chrom_regions in grouped.items():
                    chrom_data = db.load_chromosome_data(
                        sample, chrom, layers=["nucleosomes", "5mC"]
                    )
                    if chrom_data is None:
                        continue
                    for start, end, region_id in chrom_regions:
                        region_len = int(end) - int(start)
                        if region_len < 1 or region_len > max_region_bp:
                            continue
                        fiber_ids = db.get_spanning_fiber_ids(
                            chrom_data, start, end, min_coverage=min_fiber_coverage
                        )
                        if len(fiber_ids) < min_fibers:
                            continue

                        ext_start = max(0, start - flank_bp)
                        ext_end = end + flank_bp
                        ext_len = ext_end - ext_start
                        left_flank_bp = start - ext_start

                        occ_rows = []
                        met_rows = []
                        occ_ext_rows = []
                        met_ext_rows = []
                        for fiber_int_id in fiber_ids:
                            occ_rows.append(
                                _extract_nuc_occupancy(chrom_data, int(fiber_int_id), start, region_len)
                            )
                            met_rows.append(
                                _extract_5mc(
                                    chrom_data,
                                    int(fiber_int_id),
                                    start,
                                    region_len,
                                    db.methylation_threshold,
                                )
                            )
                            occ_ext_rows.append(
                                _extract_nuc_occupancy(
                                    chrom_data, int(fiber_int_id), ext_start, ext_len
                                )
                            )
                            met_ext_rows.append(
                                _extract_5mc(
                                    chrom_data,
                                    int(fiber_int_id),
                                    ext_start,
                                    ext_len,
                                    db.methylation_threshold,
                                )
                            )

                        mean_occ = np.array(occ_rows, dtype=np.float64).mean(axis=0)
                        mean_met = _nanmean_matrix(met_rows)
                        nuc_body.append(resize_vector(mean_occ, body_bins))
                        met_body.append(resize_vector(smooth_methylation(mean_met), body_bins))

                        mean_occ_ext = np.array(occ_ext_rows, dtype=np.float64).mean(axis=0)
                        mean_met_ext = _nanmean_matrix(met_ext_rows)
                        nuc_meta.append(
                            np.concatenate([
                                resize_vector(mean_occ_ext[:left_flank_bp], FLANK_BINS),
                                resize_vector(
                                    mean_occ_ext[left_flank_bp:left_flank_bp + region_len],
                                    METAPLOT_BODY_BINS,
                                ),
                                resize_vector(mean_occ_ext[left_flank_bp + region_len:], FLANK_BINS),
                            ])
                        )
                        met_meta.append(
                            np.concatenate([
                                resize_vector(mean_met_ext[:left_flank_bp], FLANK_BINS),
                                resize_vector(
                                    mean_met_ext[left_flank_bp:left_flank_bp + region_len],
                                    METAPLOT_BODY_BINS,
                                ),
                                resize_vector(mean_met_ext[left_flank_bp + region_len:], FLANK_BINS),
                            ])
                        )

                        with warnings.catch_warnings():
                            warnings.simplefilter("ignore", RuntimeWarning)
                            mean_met_value = float(np.nanmean(mean_met))
                        all_summary.append(
                            {
                                "annotation": annotation,
                                "region_id": region_id,
                                "sample": sample,
                                "chrom": chrom,
                                "start": start,
                                "end": end,
                                "region_length": region_len,
                                "n_fibers": int(len(fiber_ids)),
                                "body_bins": body_bins,
                                "mean_nuc_occupancy": float(np.mean(mean_occ)),
                                "mean_methylation": mean_met_value,
                            }
                        )

                prefix = f"{annotation}_{sample}"
                if nuc_body:
                    np.savez_compressed(
                        outdir / f"{prefix}_nuc_body.npz",
                        matrix=np.array(nuc_body, dtype=np.float64),
                        body_bins=np.array(body_bins),
                    )
                    np.savez_compressed(
                        outdir / f"{prefix}_met_body.npz",
                        matrix=np.array(met_body, dtype=np.float64),
                        body_bins=np.array(body_bins),
                    )
                    np.savez_compressed(
                        outdir / f"{prefix}_nuc_metaplot.npz",
                        matrix=np.array(nuc_meta, dtype=np.float64),
                    )
                    np.savez_compressed(
                        outdir / f"{prefix}_met_metaplot.npz",
                        matrix=np.array(met_meta, dtype=np.float64),
                    )

    _write_region_summary(outdir / "region_summary.csv", all_summary)
    return all_summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--annotations", nargs="+", required=True)
    parser.add_argument("--samples", nargs="+")
    parser.add_argument("--max-regions", type=int)
    parser.add_argument("--min-fibers", type=int, default=MIN_FIBERS)
    parser.add_argument("--min-fiber-coverage", type=float, default=MIN_FIBER_COVERAGE)
    args = parser.parse_args()

    rows = export_annotation_matrices(
        args.db,
        args.outdir,
        args.annotations,
        samples=args.samples,
        max_regions=args.max_regions,
        min_fibers=args.min_fibers,
        min_fiber_coverage=args.min_fiber_coverage,
    )
    print(f"Wrote {len(rows)} region summaries to {args.outdir}")


if __name__ == "__main__":
    main()
