"""Tests for reproducible ONT benchmark outputs."""

from __future__ import annotations

import csv
import json

from PACKAGE.benchmark.query_speed import run_query_speed_benchmark
from PACKAGE.benchmark.storage import run_storage_benchmark
from PACKAGE.db import FiberDatabase


def test_storage_benchmark_counts_database_layers(tiny_db_path_with_layers, tmp_path):
    output_dir = tmp_path / "storage"
    summary = run_storage_benchmark(tiny_db_path_with_layers, output_dir)

    assert summary["samples"]["d0"] == {
        "fibers": 2,
        "5mC": 6,
        "5hmC": 0,
        "6mA": 0,
        "msp": 3,
        "nucleosomes": 5,
    }
    assert summary["annotation_records"] == 3
    assert (output_dir / "database_record_counts.csv").stat().st_size > 0
    assert json.loads((output_dir / "storage_summary.json").read_text())["database_size_bytes"] > 0


def test_query_speed_benchmark_writes_raw_summary_and_plots(
    tiny_db_path_with_layers, tmp_path
):
    output_dir = tmp_path / "query"
    with FiberDatabase(tiny_db_path_with_layers, build_index=True) as db:
        rows = run_query_speed_benchmark(
            db,
            output_dir,
            sample="d0",
            region_sizes_kb=[1],
            n_regions=2,
            repeats=1,
            seed=7,
        )

    assert len(rows) == 4
    assert {row["mode"] for row in rows} == {"array_scan", "spatial_index"}
    with (output_dir / "query_speed_summary.csv").open() as handle:
        summary_rows = list(csv.DictReader(handle))
    assert len(summary_rows) == 2
    assert (output_dir / "query_speed_results.csv").stat().st_size > 0
    assert (output_dir / "query_speed_metadata.json").stat().st_size > 0
    assert (output_dir / "query_speed.png").stat().st_size > 0
    assert (output_dir / "query_speed.pdf").stat().st_size > 0


def test_query_speed_benchmark_can_run_spatial_index_only(
    tiny_db_path_with_layers, tmp_path
):
    output_dir = tmp_path / "query_spatial_only"
    with FiberDatabase(tiny_db_path_with_layers, build_index=True) as db:
        rows = run_query_speed_benchmark(
            db,
            output_dir,
            sample="d0",
            region_sizes_kb=[1],
            n_regions=2,
            repeats=1,
            seed=7,
            query_modes=["spatial_index"],
        )

    assert len(rows) == 2
    assert {row["mode"] for row in rows} == {"spatial_index"}
    metadata = json.loads((output_dir / "query_speed_metadata.json").read_text())
    assert metadata["requested_modes"] == ["spatial_index"]
