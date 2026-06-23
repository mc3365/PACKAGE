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
MET_BIN_SIZE = 50
CENTER_METAPLOT_BINS = 400
PROMOTER_TYPES = {"Promoter", "Bivalent_Promoter", "PRC_Promoter"}
GENEBODY_TYPES = {"Gene_Body", "Bivalent_Genebody", "PRC_Genebody"}
DIRECTIONAL_TYPES = PROMOTER_TYPES | GENEBODY_TYPES
EXTENSION_CONFIG: dict[str, int] = {
    "CGI": 2500,
    "DMV": 2500,
    "Promoter": 2500,
    "Bivalent_Promoter": 2500,
    "PRC_Promoter": 2500,
    "Typical_Enhancer": 2500,
    "Super_Enhancer": 2500,
    "FIRE_d0_Peak": 2500,
    "FIRE_d4_Peak": 2500,
    "Gene_Body": 10000,
    "Bivalent_Genebody": 10000,
    "PRC_Genebody": 10000,
}


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


def _extract_nuc_masked(
    chrom_data: dict[str, Any],
    fiber_int_id: int,
    start: int,
    length: int,
    fiber_start: int,
    fiber_end: int,
) -> np.ndarray:
    occupancy = np.full(length, np.nan, dtype=np.float64)
    cov_start = max(0, int(fiber_start) - start)
    cov_end = min(length, int(fiber_end) - start)
    if cov_start < cov_end:
        occupancy[cov_start:cov_end] = 0.0
    sl = chrom_data.get("nucleosomes_idx", {}).get(int(fiber_int_id))
    if sl is None or "nucleosomes" not in chrom_data:
        return occupancy
    s, e = sl
    data = chrom_data["nucleosomes"]
    keep = (data["ends"][s:e] > start) & (data["starts"][s:e] < start + length)
    for nuc_start, nuc_end in zip(data["starts"][s:e][keep], data["ends"][s:e][keep], strict=False):
        left = max(cov_start, int(nuc_start) - start)
        right = min(cov_end, int(nuc_end) - start)
        if right > left:
            occupancy[left:right] = 1.0
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


def _extract_5mc_masked(
    chrom_data: dict[str, Any],
    fiber_int_id: int,
    start: int,
    length: int,
    fiber_start: int,
    fiber_end: int,
    methylation_threshold: float,
) -> np.ndarray:
    methylation = np.full(length, np.nan, dtype=np.float64)
    sl = chrom_data.get("5mC_idx", {}).get(int(fiber_int_id))
    if sl is None or "5mC" not in chrom_data:
        return methylation
    s, e = sl
    data = chrom_data["5mC"]
    positions = data["positions"][s:e]
    keep = (
        (positions >= start)
        & (positions < start + length)
        & (positions >= fiber_start)
        & (positions < fiber_end)
    )
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


def _bin_methylation_fibers(
    methylation_rows: list[np.ndarray],
    window_len: int,
    bin_size: int = MET_BIN_SIZE,
) -> np.ndarray:
    n_bins = max(1, window_len // bin_size)
    binned = np.full((len(methylation_rows), n_bins), np.nan, dtype=np.float64)
    edges = np.linspace(0, window_len, n_bins + 1).astype(int)
    for row_i, values in enumerate(methylation_rows):
        for bin_i in range(n_bins):
            chunk = values[edges[bin_i]:edges[bin_i + 1]]
            if len(chunk) == 0:
                continue
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                value = np.nanmean(chunk)
            if not np.isnan(value):
                binned[row_i, bin_i] = value
    return binned


def load_strand_lookup(bed_path: Path | None) -> dict[tuple[str, int, int, str], dict[str, Any]]:
    """Load optional 9-column master BED strand metadata."""
    if bed_path is None:
        return {}
    lookup: dict[tuple[str, int, int, str], dict[str, Any]] = {}
    with Path(bed_path).open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 6:
                continue
            chrom, start, end, uid = cols[:4]
            gene_name = cols[4] if len(cols) > 4 else "."
            try:
                strand = int(cols[5])
            except ValueError:
                strand = 0
            try:
                tss = int(cols[6]) if len(cols) > 6 else -1
            except ValueError:
                tss = -1
            lookup[(chrom, int(start), int(end), uid)] = {
                "strand": strand,
                "tss": tss,
                "gene_name": gene_name,
            }
    return lookup


def _center_for_region(
    annotation: str,
    chrom: str,
    start: int,
    end: int,
    uid: str,
    strand_lookup: dict[tuple[str, int, int, str], dict[str, Any]],
) -> tuple[int, int, str]:
    info = strand_lookup.get((chrom, start, end, uid), {"strand": 0, "gene_name": "."})
    strand = int(info.get("strand", 0))
    gene_name = str(info.get("gene_name", "."))
    if annotation in PROMOTER_TYPES:
        return (end - 500 if strand == -1 else start + 500), strand, gene_name
    if annotation in GENEBODY_TYPES:
        return (start + end) // 2, strand, gene_name
    return (start + end) // 2, 0, "."


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


def export_centered_annotation_matrices(
    db_path: Path,
    outdir: Path,
    annotations: list[str],
    samples: list[str] | None = None,
    bed_path: Path | None = None,
    max_regions: int | None = None,
    min_fiber_coverage: float = MIN_FIBER_COVERAGE,
    min_fibers: int = MIN_FIBERS,
) -> list[dict[str, Any]]:
    """Export fixed-window center-based matrices matching the Part 2 heatmap logic.

    Promoters are centered at strand-aware TSS when a 9-column annotation BED is
    provided; non-promoter annotations use midpoint. Minus-strand directional regions
    are flipped before averaging so upstream is always on the left.
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    strand_lookup = load_strand_lookup(bed_path)
    all_summary: list[dict[str, Any]] = []

    with FiberDatabase(db_path) as db:
        selected_samples = samples or db.samples
        for sample in selected_samples:
            if sample not in db.samples:
                raise ValueError(f"Unknown sample {sample!r}; available: {db.samples}")

        for annotation in annotations:
            extension_bp = EXTENSION_CONFIG.get(annotation)
            if extension_bp is None:
                raise ValueError(
                    f"No centered extension configured for {annotation!r}; "
                    f"available: {sorted(EXTENSION_CONFIG)}"
                )
            window_len = 2 * extension_bp
            center_type = "TSS" if annotation in PROMOTER_TYPES else "midpoint"
            regions = db.get_annotation_regions_with_ids(annotation)
            if not regions:
                raise ValueError(f"Unknown or empty annotation {annotation!r}")
            if max_regions is not None:
                regions = regions[:max_regions]
            grouped = db.group_regions_by_chrom(regions)

            for sample in selected_samples:
                nuc_rows: list[np.ndarray] = []
                met_rows: list[np.ndarray] = []
                met_binned_rows: list[np.ndarray] = []
                region_ids: list[str] = []

                for chrom, chrom_regions in grouped.items():
                    chrom_data = db.load_chromosome_data(
                        sample, chrom, layers=["nucleosomes", "5mC"]
                    )
                    if chrom_data is None:
                        continue
                    fiber_start_map = {
                        int(fid): int(start)
                        for fid, start in zip(
                            chrom_data["fm_ids"], chrom_data["fm_starts"], strict=False
                        )
                    }
                    fiber_end_map = {
                        int(fid): int(end)
                        for fid, end in zip(
                            chrom_data["fm_ids"], chrom_data["fm_ends"], strict=False
                        )
                    }

                    for region_start, region_end, region_id in chrom_regions:
                        center, strand, gene_name = _center_for_region(
                            annotation,
                            chrom,
                            region_start,
                            region_end,
                            region_id,
                            strand_lookup,
                        )
                        window_start = center - extension_bp
                        if window_start < 0:
                            continue
                        fiber_ids = db.get_spanning_fiber_ids(
                            chrom_data,
                            region_start,
                            region_end,
                            min_coverage=min_fiber_coverage,
                        )
                        if len(fiber_ids) < min_fibers:
                            continue

                        occ_vectors: list[np.ndarray] = []
                        met_vectors: list[np.ndarray] = []
                        fiber_spans: list[int] = []
                        for fiber_int_id in fiber_ids:
                            fid = int(fiber_int_id)
                            fiber_start = fiber_start_map[fid]
                            fiber_end = fiber_end_map[fid]
                            fiber_spans.append(fiber_end - fiber_start)
                            occ = _extract_nuc_masked(
                                chrom_data, fid, window_start, window_len, fiber_start, fiber_end
                            )
                            met = _extract_5mc_masked(
                                chrom_data,
                                fid,
                                window_start,
                                window_len,
                                fiber_start,
                                fiber_end,
                                db.methylation_threshold,
                            )
                            if strand == -1 and annotation in DIRECTIONAL_TYPES:
                                occ = occ[::-1].copy()
                                met = met[::-1].copy()
                            occ_vectors.append(occ)
                            met_vectors.append(met)

                        occ_matrix = np.array(occ_vectors, dtype=np.float64)
                        met_matrix = np.array(met_vectors, dtype=np.float64)
                        with warnings.catch_warnings():
                            warnings.simplefilter("ignore", RuntimeWarning)
                            mean_occ = np.nanmean(occ_matrix, axis=0)
                            mean_met = np.nanmean(met_matrix, axis=0)
                        nuc_rows.append(mean_occ)
                        all_nan_mask = np.all(np.isnan(met_matrix), axis=0)
                        smoothed_met = smooth_methylation(mean_met)
                        smoothed_met[all_nan_mask] = np.nan
                        met_rows.append(smoothed_met)
                        binned = _bin_methylation_fibers(met_vectors, window_len)
                        with warnings.catch_warnings():
                            warnings.simplefilter("ignore", RuntimeWarning)
                            met_binned_rows.append(np.nanmean(binned, axis=0))
                            mean_met_value = float(np.nanmean(mean_met))
                            mean_occ_value = float(np.nanmean(mean_occ))
                        region_ids.append(region_id)
                        all_summary.append(
                            {
                                "annotation": annotation,
                                "region_id": region_id,
                                "sample": sample,
                                "chrom": chrom,
                                "start": region_start,
                                "end": region_end,
                                "n_fibers": int(len(fiber_ids)),
                                "region_length": int(region_end - region_start),
                                "center_type": center_type,
                                "center": int(center),
                                "extension_bp": int(extension_bp),
                                "mean_nuc_occupancy": mean_occ_value,
                                "mean_methylation": mean_met_value,
                                "n_fibers_excluded_span": 0,
                                "mean_fiber_span": float(np.mean(fiber_spans)) if fiber_spans else 0.0,
                                "strand": int(strand),
                                "gene_name": gene_name,
                            }
                        )

                if not nuc_rows:
                    continue

                nuc_matrix = np.array(nuc_rows, dtype=np.float64)
                met_matrix = np.array(met_rows, dtype=np.float64)
                met_binned_matrix = np.array(met_binned_rows, dtype=np.float64)
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", RuntimeWarning)
                    nuc_metaplot = resize_vector(np.nanmean(nuc_matrix, axis=0), CENTER_METAPLOT_BINS)
                    met_metaplot = resize_vector(np.nanmean(met_matrix, axis=0), CENTER_METAPLOT_BINS)

                prefix = f"{annotation}_{sample}"
                np.savez_compressed(
                    outdir / f"{prefix}_nuc.npz",
                    matrix=nuc_matrix,
                    region_ids=np.array(region_ids, dtype=object),
                    extension_bp=np.array(extension_bp),
                    center_type=np.array(center_type),
                )
                np.savez_compressed(
                    outdir / f"{prefix}_met.npz",
                    matrix=met_matrix,
                    region_ids=np.array(region_ids, dtype=object),
                    extension_bp=np.array(extension_bp),
                    center_type=np.array(center_type),
                )
                np.savez_compressed(
                    outdir / f"{prefix}_met_binned.npz",
                    matrix=met_binned_matrix,
                    region_ids=np.array(region_ids, dtype=object),
                    extension_bp=np.array(extension_bp),
                    center_type=np.array(center_type),
                    bin_size=np.array(MET_BIN_SIZE),
                )
                with (outdir / f"{prefix}_metaplot.csv").open("w", newline="") as handle:
                    writer = csv.DictWriter(handle, fieldnames=["bin", "nuc_mean", "met_mean"])
                    writer.writeheader()
                    writer.writerows(
                        {
                            "bin": i,
                            "nuc_mean": nuc_value,
                            "met_mean": met_value,
                        }
                        for i, (nuc_value, met_value) in enumerate(
                            zip(nuc_metaplot, met_metaplot, strict=False)
                        )
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
