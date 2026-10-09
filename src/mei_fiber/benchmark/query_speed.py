"""Reproducible regional-query performance benchmark."""

from __future__ import annotations

import csv
import json
import os
import platform
import statistics
import time
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

from mei_fiber import __version__
from mei_fiber.db import schema


def _fiber_centered_regions(
    db: Any,
    sample: str,
    region_sizes_bp: list[int],
    n_regions: int,
    seed: int,
) -> dict[int, list[tuple[str, int, int]]]:
    """Select deterministic windows centered on randomly sampled fibers."""
    chroms = db.get_chromosomes(sample)
    counts = np.array(
        [len(db.db[f"{schema.fiber_metadata_path(sample, chrom)}/starts"]) for chrom in chroms],
        dtype=np.int64,
    )
    if not len(chroms) or counts.sum() == 0:
        raise ValueError(f"Sample {sample!r} contains no fibers")

    rng = np.random.default_rng(seed)
    probabilities = counts / counts.sum()
    regions: dict[int, list[tuple[str, int, int]]] = {}
    for size in region_sizes_bp:
        selected: list[tuple[str, int, int]] = []
        for _ in range(n_regions):
            for _attempt in range(100):
                chrom = str(rng.choice(chroms, p=probabilities))
                meta = schema.fiber_metadata_path(sample, chrom)
                row = int(rng.integers(0, len(db.db[f"{meta}/starts"])))
                fiber_start = int(db.db[f"{meta}/starts"][row])
                fiber_end = int(db.db[f"{meta}/ends"][row])
                if 0 < fiber_end - fiber_start <= db.MAX_FIBER_SPAN:
                    break
            else:
                raise ValueError("Could not sample a valid fiber after 100 attempts")
            center = (fiber_start + fiber_end) // 2
            start = max(0, center - size // 2)
            selected.append((chrom, start, start + size))
        regions[size] = selected
    return regions


def _percentile(values: list[float], percentile: float) -> float:
    return float(np.percentile(np.asarray(values), percentile))


def _write_plot(summary_rows: list[dict[str, Any]], output_dir: Path) -> None:
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:  # pragma: no cover - optional dependency
        raise ImportError(
            'Benchmark plots require matplotlib. Install with `pip install ".[viz]"`.'
        ) from exc

    fig, ax = plt.subplots(figsize=(6.5, 4.2), constrained_layout=True)
    modes = sorted({str(row["mode"]) for row in summary_rows})
    colors = {"array_scan": "#4c78a8", "spatial_index": "#d65f5f"}
    for mode in modes:
        rows = sorted(
            (row for row in summary_rows if row["mode"] == mode),
            key=lambda row: int(row["region_size_bp"]),
        )
        x = np.array([int(row["region_size_bp"]) / 1000 for row in rows])
        y = np.array([float(row["median_seconds"]) * 1000 for row in rows])
        low = y - np.array([float(row["q25_seconds"]) * 1000 for row in rows])
        high = np.array([float(row["q75_seconds"]) * 1000 for row in rows]) - y
        ax.errorbar(
            x,
            y,
            yerr=np.vstack([low, high]),
            marker="o",
            linewidth=1.8,
            capsize=3,
            color=colors.get(mode),
            label=mode.replace("_", " "),
        )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Query window (kb)")
    ax.set_ylabel("Query time (ms), median and IQR")
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False)
    fig.savefig(output_dir / "query_speed.png", dpi=220)
    fig.savefig(output_dir / "query_speed.pdf")
    plt.close(fig)


def run_query_speed_benchmark(
    db: Any,
    output_dir: Path,
    sample: str | None = None,
    region_sizes_kb: list[int] | None = None,
    n_regions: int = 25,
    repeats: int = 3,
    seed: int = 3365,
    query_modes: list[str] | None = None,
) -> list[dict[str, Any]]:
    """Benchmark array-scan and available spatial-index region queries.

    Regions are centered on reproducibly sampled fibers so every timed query overlaps
    real data. Each selected region is queried ``repeats`` times. By default, the
    function benchmarks the array scan and adds the spatial-index mode when an index is
    loaded. Large production databases can request only ``spatial_index`` to avoid
    hundreds of intentionally slow full-array scans.

    Outputs:
        ``query_speed_results.csv``
            One row per timed query.
        ``query_speed_summary.csv``
            Median, quartiles, and p95 for each mode and region size.
        ``query_speed_metadata.json``
            Reproducibility metadata for the run.
        ``query_speed.png`` and ``query_speed.pdf``
            Plot of median query time with interquartile ranges.
    """
    if n_regions < 1 or repeats < 1:
        raise ValueError("n_regions and repeats must both be at least 1")
    region_sizes_kb = region_sizes_kb or [1, 10, 100, 1000]
    if any(size < 1 for size in region_sizes_kb):
        raise ValueError("region sizes must be positive integers in kb")
    query_modes = query_modes or ["array_scan", "spatial_index"]
    allowed_modes = {"array_scan", "spatial_index"}
    unknown_modes = sorted(set(query_modes) - allowed_modes)
    if unknown_modes:
        raise ValueError(f"Unknown query mode(s): {unknown_modes}; expected {sorted(allowed_modes)}")
    sample = sample or db.samples[0]
    if sample not in db.samples:
        raise ValueError(f"Unknown sample {sample!r}; available: {db.samples}")

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    sizes_bp = [size * 1000 for size in region_sizes_kb]
    regions = _fiber_centered_regions(db, sample, sizes_bp, n_regions, seed)

    original_index = db._spatial_index
    modes: list[tuple[str, dict[str, Any]]] = []
    if "array_scan" in query_modes:
        modes.append(("array_scan", {}))
    if "spatial_index" in query_modes and sample in original_index and original_index[sample]:
        modes.append(("spatial_index", original_index))
    if "spatial_index" in query_modes and not any(mode == "spatial_index" for mode, _ in modes):
        raise ValueError(
            "Spatial-index mode requested, but no spatial index is loaded for "
            f"sample {sample!r}. Re-run with --build-index or omit spatial_index."
        )
    if not modes:
        raise ValueError("No query modes are available to benchmark")

    result_rows: list[dict[str, Any]] = []
    reference_counts: dict[tuple[int, int], int] = {}
    try:
        for mode, index in modes:
            print(f"Benchmarking query mode: {mode}", flush=True)
            db._spatial_index = index
            # Warm every selected region once. OS-level cold-cache control requires
            # administrative privileges and is not portable across HPC systems.
            for size in sizes_bp:
                print(f"  Warming {size // 1000} kb regions", flush=True)
                for region_number, (chrom, start, end) in enumerate(regions[size], 1):
                    count = len(db.get_fibers_at(chrom, start, end, sample=sample))
                    key = (size, region_number)
                    if mode == "array_scan":
                        reference_counts[key] = count
                    elif key in reference_counts and count != reference_counts[key]:
                        raise RuntimeError(
                            "Spatial index and array scan disagree for "
                            f"{chrom}:{start}-{end}: {count} != {reference_counts[key]}. "
                            "Rebuild the spatial index before benchmarking."
                        )
            for size in sizes_bp:
                print(f"  Timing {size // 1000} kb regions", flush=True)
                for region_number, (chrom, start, end) in enumerate(regions[size], 1):
                    for repeat in range(1, repeats + 1):
                        began = time.perf_counter()
                        fibers = db.get_fibers_at(chrom, start, end, sample=sample)
                        elapsed = time.perf_counter() - began
                        result_rows.append(
                            {
                                "sample": sample,
                                "mode": mode,
                                "region_size_bp": size,
                                "region_number": region_number,
                                "repeat": repeat,
                                "chrom": chrom,
                                "start": start,
                                "end": end,
                                "n_fibers": len(fibers),
                                "elapsed_seconds": f"{elapsed:.9f}",
                            }
                        )
    finally:
        db._spatial_index = original_index

    result_path = output_dir / "query_speed_results.csv"
    with result_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(result_rows[0]))
        writer.writeheader()
        writer.writerows(result_rows)

    grouped: dict[tuple[str, int], list[float]] = defaultdict(list)
    for row in result_rows:
        grouped[(str(row["mode"]), int(row["region_size_bp"]))].append(
            float(row["elapsed_seconds"])
        )
    summary_rows: list[dict[str, Any]] = []
    for (mode, size), values in sorted(grouped.items()):
        summary_rows.append(
            {
                "sample": sample,
                "mode": mode,
                "region_size_bp": size,
                "n_queries": len(values),
                "median_seconds": f"{statistics.median(values):.9f}",
                "q25_seconds": f"{_percentile(values, 25):.9f}",
                "q75_seconds": f"{_percentile(values, 75):.9f}",
                "p95_seconds": f"{_percentile(values, 95):.9f}",
            }
        )
    with (output_dir / "query_speed_summary.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary_rows[0]))
        writer.writeheader()
        writer.writerows(summary_rows)

    metadata = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "mei_fiber_version": __version__,
        "database": str(db.db_path.resolve()),
        "database_size_bytes": db.db_path.stat().st_size,
        "sample": sample,
        "modes": [mode for mode, _ in modes],
        "requested_modes": query_modes,
        "region_sizes_kb": region_sizes_kb,
        "n_regions_per_size": n_regions,
        "repeats_per_region": repeats,
        "seed": seed,
        "cache_condition": "warm",
        "python": platform.python_version(),
        "platform": platform.platform(),
        "processor": platform.processor(),
        "logical_cpus": os.cpu_count(),
    }
    (output_dir / "query_speed_metadata.json").write_text(
        json.dumps(metadata, indent=2) + "\n"
    )
    _write_plot(summary_rows, output_dir)
    return result_rows
