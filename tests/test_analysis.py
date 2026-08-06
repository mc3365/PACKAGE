"""Tests for analysis exports used by visualization examples."""

from __future__ import annotations

import csv
import json

import numpy as np

from PACKAGE.analysis import (
    cov_to_object,
    export_annotation_matrices,
    export_centered_annotation_matrices,
    export_coaccessibility_cov,
    export_global_feature_fractions,
    find_contained_peak_pairs,
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


def test_smooth_methylation_does_not_zero_pad_edges():
    values = np.ones(500, dtype=float)

    smoothed = smooth_methylation(values, sigma=50)

    assert np.allclose(smoothed[[0, -1]], 1.0)


def test_find_contained_peak_pairs_matches_bedtools_f1_logic(tmp_path):
    stitched = tmp_path / "FIRE_stitched.bed"
    peaks = tmp_path / "FIRE_peaks_intergenic.bed"
    stitched.write_text("chr1\t100\t800\nchr1\t1000\t1400\n")
    peaks.write_text("chr1\t120\t180\nchr1\t200\t900\nchr1\t300\t350\nchr1\t650\t760\nchr1\t900\t950\n")

    pairs = find_contained_peak_pairs(stitched, peaks)

    assert [(p.element.region_id, p.stitched.region_id) for p in pairs] == [
        ("chr1:120-180", "chr1:100-800"),
        ("chr1:300-350", "chr1:100-800"),
        ("chr1:650-760", "chr1:100-800"),
    ]


def test_export_coaccessibility_cov_writes_legacy_shape(tiny_db_path_with_layers, tmp_path):
    stitched = tmp_path / "FIRE_stitched.bed"
    peaks = tmp_path / "FIRE_peaks_intergenic.bed"
    out = tmp_path / "Cov.bed"
    stitched.write_text("chr1\t100\t850\n")
    peaks.write_text("chr1\t120\t180\nchr1\t650\t760\n")

    summary = export_coaccessibility_cov(
        tiny_db_path_with_layers,
        stitched,
        peaks,
        out,
        sample="d0",
    )

    rows = [line.split("\t") for line in out.read_text().splitlines()]
    assert summary["contained_peak_pairs"] == 2
    assert summary["cov_rows"] == 2
    assert rows == [
        ["chr1", "120", "180", "chr1", "100", "850", "fiber1", "0.05", "60"],
        ["chr1", "650", "760", "chr1", "100", "850", "fiber1", "1.01", "110"],
    ]


def test_cov_to_object_matches_legacy_obj_shape(tmp_path):
    cov = tmp_path / "Cov.bed"
    out = tmp_path / "scored_obj.json"
    cov.write_text(
        "chr1\t120\t180\tchr1\t100\t850\tfiber1\t0.05\t60\n"
        "chr1\t650\t760\tchr1\t100\t850\tfiber1\t1.01\t110\n"
        "chr1\t650\t760\tchr1\t100\t850\tfiber2\t0.25\t20\n"
    )

    summary = cov_to_object(cov, out)
    data = json.loads(out.read_text())

    assert summary["stitched_regions"] == 1
    assert data[0]["seId"] == "chr1:100-850"
    assert data[0]["enhs"][0]["enhId"] == "chr1:120-180"
    assert data[0]["enhs"][0]["fibers"] == [0.05, None]
    assert data[0]["enhs"][1]["enhId"] == "chr1:650-760"
    assert data[0]["enhs"][1]["fibers"] == [1.01, None]
