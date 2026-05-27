"""Tests for the database layer.

These start small. As we refactor the FiberDatabase methods in, we add a test for each
behavior (open, list samples, list chromosomes, query by region, ...).
"""

from __future__ import annotations

import pytest

from PACKAGE.db import FiberDatabase


def test_database_opens(tiny_db_path):
    """FiberDatabase should open a well-formed file and report its samples."""
    with FiberDatabase(tiny_db_path) as db:
        assert db.samples == ["d0"]
        assert db.version == "8.0"
        assert db.genome_version == "mm10"
        assert db.methylation_threshold == 0.5


def test_database_missing_file_raises(tmp_path):
    """FiberDatabase should raise FileNotFoundError for a missing file."""
    with pytest.raises(FileNotFoundError):
        FiberDatabase(tmp_path / "does_not_exist.h5")


def test_get_chromosomes(tiny_db_path):
    """get_chromosomes should return chr-prefixed groups only."""
    with FiberDatabase(tiny_db_path) as db:
        assert db.get_chromosomes("d0") == ["chr1"]


def test_summary(tiny_db_path):
    """get_summary should return a populated dict."""
    with FiberDatabase(tiny_db_path) as db:
        summary = db.get_summary()
    assert summary["samples"] == ["d0"]
    assert summary["genome"] == "mm10"
    assert summary["version"] == "8.0"
