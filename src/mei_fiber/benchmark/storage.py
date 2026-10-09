"""Database scale and storage benchmark."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

import h5py

from mei_fiber.config import load_config
from mei_fiber.db import schema


def run_storage_benchmark(
    db_path: Path,
    output_dir: Path,
    config_path: Path | None = None,
) -> dict[str, Any]:
    """Count stored records and compare HDF5 size with configured intermediates."""
    db_path = Path(db_path)
    if not db_path.exists():
        raise FileNotFoundError(f"Database not found: {db_path}")
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    count_rows: list[dict[str, Any]] = []
    sample_totals: dict[str, dict[str, int]] = {}
    with h5py.File(db_path, "r") as db:
        samples = [name for name in db if name not in schema.RESERVED_TOP_LEVEL]
        for sample in samples:
            totals = {"fibers": 0, **dict.fromkeys(schema.SUPPORTED_LAYERS, 0)}
            chroms = sorted(name for name in db[sample] if name.startswith("chr"))
            for chrom in chroms:
                meta = schema.fiber_metadata_path(sample, chrom)
                n_fibers = len(db[f"{meta}/fiber_int_ids"])
                totals["fibers"] += n_fibers
                count_rows.append(
                    {"sample": sample, "chrom": chrom, "layer": "fibers", "records": n_fibers}
                )
                for layer in schema.SUPPORTED_LAYERS:
                    path = schema.layer_path(sample, chrom, layer)
                    if path not in db:
                        continue
                    n_records = len(db[f"{path}/fiber_int_ids"])
                    totals[layer] += n_records
                    count_rows.append(
                        {
                            "sample": sample,
                            "chrom": chrom,
                            "layer": layer,
                            "records": n_records,
                        }
                    )
            sample_totals[sample] = totals

        annotation_count = 0
        annotation_root = schema.annotation_features_path()
        if annotation_root in db:
            annotation_count = sum(len(dataset) for dataset in db[annotation_root].values())

    with (output_dir / "database_record_counts.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["sample", "chrom", "layer", "records"])
        writer.writeheader()
        writer.writerows(count_rows)

    input_files: list[dict[str, Any]] = []
    if config_path is not None:
        config = load_config(config_path)
        paths = {path for sample in config.samples for path in sample.layers.values()}
        for path in sorted(paths):
            input_files.append(
                {
                    "path": str(path),
                    "exists": path.exists(),
                    "size_bytes": path.stat().st_size if path.exists() else None,
                }
            )
        with (output_dir / "intermediate_file_sizes.csv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=["path", "exists", "size_bytes"])
            writer.writeheader()
            writer.writerows(input_files)

    index_path = db_path.with_suffix(".index.pkl")
    summary = {
        "database": str(db_path.resolve()),
        "database_size_bytes": db_path.stat().st_size,
        "spatial_index_size_bytes": index_path.stat().st_size if index_path.exists() else 0,
        "annotation_records": annotation_count,
        "samples": sample_totals,
        "intermediate_size_bytes": sum(
            int(row["size_bytes"]) for row in input_files if row["size_bytes"] is not None
        ),
        "missing_intermediate_files": [row["path"] for row in input_files if not row["exists"]],
    }
    (output_dir / "storage_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    return summary

