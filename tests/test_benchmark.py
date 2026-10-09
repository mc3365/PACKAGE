"""Tests for reproducible ONT benchmark outputs."""

from __future__ import annotations

import csv
import json

from PACKAGE.benchmark.annotation import run_annotation_benchmark
from PACKAGE.benchmark.coaccessibility import run_coaccessibility_validation
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
        "fire_accessibility": 3,
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


def test_annotation_benchmark_times_real_annotation_query(
    tiny_db_path_with_layers, tmp_path
):
    output_dir = tmp_path / "annotation"
    rows = run_annotation_benchmark(
        tiny_db_path_with_layers,
        output_dir,
        sample="d0",
        annotations=["Typical_Enhancer"],
        layers=["nucleosomes", "5mC"],
        max_regions=1,
        repeats=1,
        warmups=0,
    )

    assert len(rows) == 1
    assert rows[0]["annotation"] == "Typical_Enhancer"
    assert int(rows[0]["n_rows"]) > 0
    assert (output_dir / "annotation_query_results.csv").stat().st_size > 0
    assert (output_dir / "annotation_query_summary.csv").stat().st_size > 0
    metadata = json.loads((output_dir / "annotation_query_metadata.json").read_text())
    assert metadata["query_function"] == "FiberDatabase.query_annotation_fast"


def test_coaccessibility_validation_normalizes_directional_pairs(tmp_path):
    legacy_ce = tmp_path / "legacy_ce.txt"
    package_ce = tmp_path / "package_ce.txt"
    legacy_cluster = tmp_path / "legacy_cluster.txt"
    package_cluster = tmp_path / "package_cluster.txt"

    legacy_ce.write_text(
        "chr1:10-20\tchr1:30-40\tchr1:0-100\t4.0\t1\tSuper\n"
        "chr1:30-40\tchr1:10-20\tchr1:0-100\t4.0\t2\tSuper\n"
        "chr1:110-120\tchr1:130-140\tchr1:100-200\t2.0\t3\t\n"
        "chr1:130-140\tchr1:110-120\tchr1:100-200\t2.0\t4\t\n"
    )
    package_ce.write_text(
        "chr1:10-20\tchr1:30-40\tchr1:0-100\t5.0\t1\tSuper\n"
        "chr1:110-120\tchr1:130-140\tchr1:100-200\t1.5\t2\t\n"
    )
    legacy_cluster.write_text(
        "chr1:10-20, chr1:30-40\tchr1:0-100\t2\t4.0\t1\tSuper\t3\t1\t1\t2\n"
    )
    package_cluster.write_text(
        "chr1:10-20, chr1:30-40\tchr1:0-100\t2\t5.0\t1\tSuper\t3\t1\t1\t2\n"
        "chr1:110-120, chr1:130-140\tchr1:100-200\t2\t1.5\t2\t\t1\t1\t1\t3\n"
    )

    legacy_cov = tmp_path / "legacy_cov.bed"
    package_cov = tmp_path / "package_cov.bed"
    difference = tmp_path / "difference.bed"
    package_cov.write_text("row-b\nrow-a\n")
    difference.write_text("row-c\n")
    legacy_cov.write_text("row-c\nrow-a\nrow-b\n")

    output_dir = tmp_path / "validation"
    summary = run_coaccessibility_validation(
        legacy_ce,
        package_ce,
        legacy_cluster,
        package_cluster,
        output_dir,
        legacy_cov=legacy_cov,
        package_cov=package_cov,
        cov_difference=difference,
    )

    assert summary["pairs"]["legacy_rows"] == 4
    assert summary["pairs"]["legacy_unique_unordered_pairs"] == 2
    assert summary["pairs"]["legacy_directional_duplicate_rows"] == 2
    assert summary["pairs"]["shared_pairs"] == 2
    assert summary["clusters"]["package_only_clusters"] == 1
    assert summary["clusters"]["package_only_half_jaccard_sentinel_candidates"] == 1
    assert summary["cov"]["package_plus_difference_matches_legacy"] is True
    assert (output_dir / "validation_summary.json").exists()
    assert (output_dir / "pair_rank_agreement.csv").exists()
    assert (output_dir / "cluster_rank_agreement.csv").exists()
    assert (output_dir / "cluster_membership.csv").exists()
    assert (output_dir / "coaccess_validation.png").exists()
