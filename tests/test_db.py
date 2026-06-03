"""Tests for the database layer (read-only query API).

These tests use the small synthetic HDF5 from ``conftest.tiny_db_path``. As Batch 2
and Batch 3 of the port land, we add per-method tests here. Right now (post-Batch-1)
we test: open/close, metadata reading, sample/chromosome enumeration, fiber-count
summary, and the backward-compat annotation path resolution.
"""

from __future__ import annotations

from pathlib import Path

import h5py
import pytest

from PACKAGE.db import FiberDatabase


# ---------------------------------------------------------------------------
# Basic open/close + metadata
# ---------------------------------------------------------------------------
def test_database_opens(tiny_db_path):
    """FiberDatabase should open a well-formed file and expose top-level metadata."""
    with FiberDatabase(tiny_db_path) as db:
        assert db.samples == ["d0"]
        assert db.version == "8.0"
        assert db.genome_version == "mm10"
        assert db.methylation_threshold == 0.5


def test_database_missing_file_raises(tmp_path):
    """FiberDatabase should raise FileNotFoundError for a non-existent path."""
    with pytest.raises(FileNotFoundError):
        FiberDatabase(tmp_path / "does_not_exist.h5")


def test_database_close_is_idempotent_via_context_manager(tiny_db_path):
    """Using the context manager twice on the same file should not crash."""
    with FiberDatabase(tiny_db_path) as db:
        _ = db.samples
    # Reopen in a second context — must work cleanly.
    with FiberDatabase(tiny_db_path) as db:
        _ = db.samples


# ---------------------------------------------------------------------------
# Sample / chromosome enumeration
# ---------------------------------------------------------------------------
def test_get_chromosomes(tiny_db_path):
    """get_chromosomes should return only groups starting with 'chr'."""
    with FiberDatabase(tiny_db_path) as db:
        assert db.get_chromosomes("d0") == ["chr1"]


def test_samples_excludes_reserved_top_level(tmp_path):
    """Top-level groups in RESERVED_TOP_LEVEL must not be listed as samples."""
    db_path = tmp_path / "withreserved.h5"
    with h5py.File(db_path, "w") as f:
        f.create_group("metadata")
        f.create_group("annotations")
        f.create_group("realsample")
    with FiberDatabase(db_path) as db:
        assert db.samples == ["realsample"]


# ---------------------------------------------------------------------------
# Summary + fiber counting
# ---------------------------------------------------------------------------
def test_summary(tiny_db_path):
    """get_summary should return populated metadata and per-sample fiber counts."""
    with FiberDatabase(tiny_db_path) as db:
        summary = db.get_summary()
    assert summary["samples"] == ["d0"]
    assert summary["genome"] == "mm10"
    assert summary["version"] == "8.0"
    # The tiny fixture has 2 fibers under d0/fiber_lookup/fiber_ids.
    assert summary["n_fibers_d0"] == 2
    # Spatial index is not built for the tiny fixture.
    assert summary["has_spatial_index"] is False


# ---------------------------------------------------------------------------
# Backward compatibility with V8 (`master_v3/`) annotation path
# ---------------------------------------------------------------------------
def _build_db_with_annotations(tmp_path: Path, anno_root: str) -> Path:
    """Build a minimal db with annotations at the given HDF5 path root."""
    db_path = tmp_path / f"db_{anno_root.replace('/', '_')}.h5"
    with h5py.File(db_path, "w") as f:
        meta = f.create_group("metadata")
        meta.attrs["version"] = "8.0"
        f.create_group("d0")
        anno_group = f.create_group(f"{anno_root}/Typical_Enhancer")
        # Make it look like a real feature table, even if empty.
        anno_group.attrs["dummy"] = True
    return db_path


def test_list_annotations_new_schema_path(tmp_path):
    """list_annotations finds features under the new 'annotations/master/...' path."""
    db_path = _build_db_with_annotations(tmp_path, "annotations/master/features")
    with FiberDatabase(db_path) as db:
        assert "Typical_Enhancer" in db.list_annotations()


def test_list_annotations_v8_legacy_path(tmp_path):
    """list_annotations falls back to the V8 'annotations/master_v3/...' path."""
    db_path = _build_db_with_annotations(tmp_path, "annotations/master_v3/features")
    with FiberDatabase(db_path) as db:
        assert "Typical_Enhancer" in db.list_annotations()


def test_list_annotations_returns_empty_when_no_annotations(tmp_path):
    """list_annotations returns [] if neither path exists."""
    db_path = tmp_path / "no_anno.h5"
    with h5py.File(db_path, "w") as f:
        f.create_group("d0")
    with FiberDatabase(db_path) as db:
        assert db.list_annotations() == []


# ---------------------------------------------------------------------------
# Internal fiber ID mapping (covered indirectly by Batch 2/3, but worth a smoke test)
# ---------------------------------------------------------------------------
def test_internal_id_map_roundtrip(tiny_db_path):
    """Round-tripping a fiber ID through str->int->str should give the same string."""
    with FiberDatabase(tiny_db_path) as db:
        int_id = db._str_to_int("fiber1", "d0", "chr1")
        assert int_id is not None
        assert db._int_to_str(int_id, "d0", "chr1") == "fiber1"


def test_internal_id_map_missing_fiber_returns_none(tiny_db_path):
    """_str_to_int returns None for a fiber not present on the chromosome."""
    with FiberDatabase(tiny_db_path) as db:
        assert db._str_to_int("not_a_real_fiber", "d0", "chr1") is None
