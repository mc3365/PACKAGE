"""Tests for analysis exports used by visualization examples."""

from __future__ import annotations

import csv

import numpy as np

from PACKAGE.analysis import (
    export_annotation_matrices,
    export_centered_annotation_matrices,
    export_global_feature_fractions,
)
from PACKAGE.analysis.heatmap import smooth_methylation


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


def test_export_centered_annotation_matrices_writes_part2_style_outputs(
    tiny_db_path_with_layers, tmp_path
):
    summary = export_centered_annotation_matrices(
        tiny_db_path_with_layers,
        tmp_path,
        annotations=["Typical_Enhancer"],
        samples=["d0"],
        min_fibers=1,
    )

    assert len(summary) == 1
    nuc = np.load(tmp_path / "Typical_Enhancer_d0_nuc.npz", allow_pickle=True)
    met = np.load(tmp_path / "Typical_Enhancer_d0_met.npz", allow_pickle=True)
    binned = np.load(tmp_path / "Typical_Enhancer_d0_met_binned.npz", allow_pickle=True)
    assert nuc["matrix"].shape == (1, 5000)
    assert met["matrix"].shape == (1, 5000)
    assert binned["matrix"].shape == (1, 100)
    assert int(nuc["extension_bp"]) == 2500
    assert (tmp_path / "Typical_Enhancer_d0_metaplot.csv").stat().st_size > 0


def test_smooth_methylation_preserves_short_vector_length():
    values = np.full(219, np.nan)
    values[[10, 40, 100, 180]] = [0, 1, 0, 1]

    smoothed = smooth_methylation(values, sigma=50)

    assert len(smoothed) == len(values)
    assert np.isnan(smoothed[np.isnan(values)]).all()
