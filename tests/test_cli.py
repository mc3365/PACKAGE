"""Tests for the CLI entry points."""

from __future__ import annotations

import json

from click.testing import CliRunner

from mei_fiber.cli import cli


def test_cli_help():
    """Top-level --help should succeed and mention MEI-Fiber."""
    runner = CliRunner()
    result = runner.invoke(cli, ["--help"])
    assert result.exit_code == 0
    assert "mei-fiber" in result.output.lower()


def test_cli_version():
    """--version should print a semver-ish string."""
    runner = CliRunner()
    result = runner.invoke(cli, ["--version"])
    assert result.exit_code == 0
    # Format is "mei_fiber, version X.Y.Z"
    assert "version" in result.output.lower()


def test_info_command(tiny_db_path):
    """`mei_fiber info <path>` should run and report samples."""
    runner = CliRunner()
    result = runner.invoke(cli, ["info", str(tiny_db_path)])
    assert result.exit_code == 0
    assert "d0" in result.output


def test_build_command_can_build_spatial_index_from_config(tmp_path, monkeypatch):
    """`mei_fiber build` should honor build.build_spatial_index from the YAML."""
    cfg_path = tmp_path / "config.yaml"
    cfg_path.write_text(
        f"""
output_dir: {tmp_path}
output_file: fibers.h5
samples:
  - name: d0
    layers:
      nucleosomes: {tmp_path / "nuc.csv"}
build:
  build_spatial_index: true
"""
    )
    called = {}

    def fake_build_database(cfg, samples):
        called["build_samples"] = samples
        called["output_path"] = cfg.output_path

    class FakeFiberDatabase:
        def __init__(self, db_path, build_index=False):
            called["db_path"] = db_path
            called["build_index"] = build_index

        def __enter__(self):
            return self

        def __exit__(self, *exc):
            return None

    monkeypatch.setattr("mei_fiber.db.builder.build_database", fake_build_database)
    monkeypatch.setattr("mei_fiber.db.database.FiberDatabase", FakeFiberDatabase)

    runner = CliRunner()
    result = runner.invoke(cli, ["build", "--config", str(cfg_path)])

    assert result.exit_code == 0, result.output
    assert called["build_samples"] == ["d0"]
    assert called["db_path"] == tmp_path / "fibers.h5"
    assert called["build_index"] is True


def test_coaccess_prepare_command(tmp_path):
    peaks = tmp_path / "peaks.bed"
    genes = tmp_path / "genes.gff3"
    chrom_sizes = tmp_path / "chrom.sizes"
    outdir = tmp_path / "prepared"
    peaks.write_text("chr1\t100\t200\nchr1\t250\t300\n")
    genes.write_text("##gff-version 3\n")
    chrom_sizes.write_text("chr1\t1000\n")

    runner = CliRunner()
    result = runner.invoke(
        cli,
        [
            "coaccess",
            "prepare",
            "--peaks",
            str(peaks),
            "--genes",
            str(genes),
            "--chrom-sizes",
            str(chrom_sizes),
            "--outdir",
            str(outdir),
            "--stitch-distance",
            "100",
        ],
    )

    assert result.exit_code == 0, result.output
    summary = json.loads(result.output)
    assert summary["counts"]["stitched_regions"] == 1
    assert (outdir / "FIRE_stitched.bed").read_text() == "chr1\t100\t300\n"
