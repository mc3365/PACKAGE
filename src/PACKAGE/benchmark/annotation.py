"""Benchmark biologically meaningful annotation queries."""

from __future__ import annotations

import argparse
import csv
import gc
import json
import platform
import statistics
import time
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

from PACKAGE import __version__
from PACKAGE.db import FiberDatabase


def _percentile(values: list[float], percentile: float) -> float:
    return float(np.percentile(np.asarray(values), percentile))


def run_annotation_benchmark(
    db_path: Path,
    output_dir: Path,
    sample: str,
    annotations: list[str],
    layers: list[str] | None = None,
    max_regions: int | None = None,
    repeats: int = 3,
    warmups: int = 1,
) -> list[dict[str, Any]]:
    """Time ``FiberDatabase.query_annotation_fast`` for real annotation classes.

    This benchmark complements the random-region benchmark by exercising the workflow
    used in biological analyses: annotation class -> overlapping fibers -> per-fiber
    layer summaries. It records timing and row counts only; it does not write the full
    query result tables, which can be large.
    """
    if repeats < 1:
        raise ValueError("repeats must be at least 1")
    if warmups < 0:
        raise ValueError("warmups must be non-negative")
    if not annotations:
        raise ValueError("at least one annotation is required")

    db_path = Path(db_path)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    layers = layers or ["nucleosomes", "5mC"]

    rows: list[dict[str, Any]] = []
    with FiberDatabase(db_path) as db:
        if sample not in db.samples:
            raise ValueError(f"Unknown sample {sample!r}; available: {db.samples}")
        available_annotations = db.list_annotations()

        for annotation in annotations:
            if annotation not in available_annotations:
                raise ValueError(
                    f"Unknown annotation {annotation!r}; available examples: "
                    f"{available_annotations[:10]}"
                )

            for warmup in range(1, warmups + 1):
                print(f"Warming {annotation} annotation query ({warmup}/{warmups})", flush=True)
                df = db.query_annotation_fast(
                    annotation,
                    sample=sample,
                    feature_types=layers,
                    max_regions=max_regions,
                )
                del df
                gc.collect()

            for repeat in range(1, repeats + 1):
                print(f"Timing {annotation} annotation query ({repeat}/{repeats})", flush=True)
                began = time.perf_counter()
                df = db.query_annotation_fast(
                    annotation,
                    sample=sample,
                    feature_types=layers,
                    max_regions=max_regions,
                )
                elapsed = time.perf_counter() - began
                rows.append(
                    {
                        "sample": sample,
                        "annotation": annotation,
                        "layers": ",".join(layers),
                        "max_regions": max_regions if max_regions is not None else "",
                        "repeat": repeat,
                        "n_rows": len(df),
                        "n_columns": len(df.columns),
                        "elapsed_seconds": f"{elapsed:.9f}",
                    }
                )
                del df
                gc.collect()

    result_path = output_dir / "annotation_query_results.csv"
    with result_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    grouped: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        grouped[(str(row["sample"]), str(row["annotation"]))].append(row)

    summary_rows: list[dict[str, Any]] = []
    for (row_sample, annotation), group_rows in sorted(grouped.items()):
        values = [float(row["elapsed_seconds"]) for row in group_rows]
        n_rows_values = [int(row["n_rows"]) for row in group_rows]
        summary_rows.append(
            {
                "sample": row_sample,
                "annotation": annotation,
                "layers": group_rows[0]["layers"],
                "max_regions": group_rows[0]["max_regions"],
                "repeats": len(group_rows),
                "median_seconds": f"{statistics.median(values):.9f}",
                "q25_seconds": f"{_percentile(values, 25):.9f}",
                "q75_seconds": f"{_percentile(values, 75):.9f}",
                "p95_seconds": f"{_percentile(values, 95):.9f}",
                "n_rows_min": min(n_rows_values),
                "n_rows_max": max(n_rows_values),
            }
        )
    with (output_dir / "annotation_query_summary.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary_rows[0]))
        writer.writeheader()
        writer.writerows(summary_rows)

    metadata = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "package_version": __version__,
        "database": str(db_path.resolve()),
        "database_size_bytes": db_path.stat().st_size,
        "sample": sample,
        "annotations": annotations,
        "layers": layers,
        "max_regions": max_regions,
        "repeats": repeats,
        "warmups": warmups,
        "query_function": "FiberDatabase.query_annotation_fast",
        "cache_condition": "warm" if warmups else "unwarmed",
        "python": platform.python_version(),
        "platform": platform.platform(),
        "processor": platform.processor(),
    }
    (output_dir / "annotation_query_metadata.json").write_text(
        json.dumps(metadata, indent=2) + "\n"
    )
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db", type=Path, required=True, help="PACKAGE HDF5 database")
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument(
        "--annotation",
        action="append",
        required=True,
        help="annotation class to benchmark; repeat for multiple classes",
    )
    parser.add_argument(
        "--layers",
        nargs="+",
        default=["nucleosomes", "5mC"],
        help="layers passed to query_annotation_fast",
    )
    parser.add_argument("--max-regions", type=int)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--warmups", type=int, default=1)
    args = parser.parse_args()

    run_annotation_benchmark(
        args.db,
        args.outdir,
        sample=args.sample,
        annotations=args.annotation,
        layers=args.layers,
        max_regions=args.max_regions,
        repeats=args.repeats,
        warmups=args.warmups,
    )
    print(f"Wrote annotation benchmark outputs to {args.outdir}")


if __name__ == "__main__":
    main()
