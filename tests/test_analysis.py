"""Tests for analysis exports used by visualization examples."""

from __future__ import annotations

import csv

import numpy as np

from PACKAGE.analysis import export_annotation_matrices, export_global_feature_fractions


def test_export_global_feature_fractions_writes_expected_columns(
    tiny_db_path_with_layers, tmp_path
):
    out_csv = tmp_path / "features.csv"
    export_global_feature_fractions(tiny_db_path_with_layers, out_csv)

    with out_csv.open() as handle:
        rows = list(csv.DictReader(handle))

    assert len(rows) == 2
    assert rows[0]["sample"] == "d0"
    assert float(rows[0]["frac_nuc_coverage"]) > 0
    assert float(rows[0]["frac_5mC_among_cpg"]) == 0.5


def test_export_annotation_matrices_writes_npz_outputs(tiny_db_path_with_layers, tmp_path):
    summary = export_annotation_matrices(
        tiny_db_path_with_layers,
        tmp_path,
        annotations=["CGI"],
        samples=["d0"],
        min_fibers=1,
    )

    assert len(summary) == 1
    nuc_body = np.load(tmp_path / "CGI_d0_nuc_body.npz")["matrix"]
    met_body = np.load(tmp_path / "CGI_d0_met_body.npz")["matrix"]
    nuc_meta = np.load(tmp_path / "CGI_d0_nuc_metaplot.npz")["matrix"]
    assert nuc_body.shape == (1, 1000)
    assert met_body.shape == (1, 1000)
    assert nuc_meta.shape == (1, 400)
    assert (tmp_path / "region_summary.csv").stat().st_size > 0
