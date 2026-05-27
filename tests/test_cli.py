"""Tests for the CLI entry points."""

from __future__ import annotations

from click.testing import CliRunner

from PACKAGE.cli import cli


def test_cli_help():
    """Top-level --help should succeed and mention 'PACKAGE'."""
    runner = CliRunner()
    result = runner.invoke(cli, ["--help"])
    assert result.exit_code == 0
    assert "package" in result.output.lower()


def test_cli_version():
    """--version should print a semver-ish string."""
    runner = CliRunner()
    result = runner.invoke(cli, ["--version"])
    assert result.exit_code == 0
    # Format is "PACKAGE, version X.Y.Z"
    assert "version" in result.output.lower()


def test_info_command(tiny_db_path):
    """`PACKAGE info <path>` should run and report samples."""
    runner = CliRunner()
    result = runner.invoke(cli, ["info", str(tiny_db_path)])
    assert result.exit_code == 0
    assert "d0" in result.output
