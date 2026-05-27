"""Shared pytest fixtures.

Fixtures defined here are auto-discovered by pytest — no import needed in test files.

Currently provides a tiny synthetic HDF5 database for testing the query layer without
needing real Fiber-seq data on disk. When the refactor lands, expand this fixture to
include realistic minimum data so analysis tests can run end-to-end.
"""

from __future__ import annotations

from pathlib import Path

import h5py
import numpy as np
import pytest


@pytest.fixture
def tiny_db_path(tmp_path: Path) -> Path:
    """Build a minimal PACKAGE-shaped HDF5 file in a per-test tmpdir.

    Currently contains just the top-level structure + metadata, enough to verify
    FiberDatabase opens correctly. Expand as the refactor adds more methods that
    need testing.
    """
    db_path = tmp_path / "tiny.h5"
    with h5py.File(db_path, "w") as f:
        # Metadata
        meta = f.create_group("metadata")
        meta.attrs["version"] = "8.0"
        meta.attrs["genome_version"] = "mm10"
        meta.attrs["methylation_threshold"] = 0.5

        # One sample with one chromosome and a few fibers (fake but schema-compliant)
        sample = f.create_group("d0")

        # fiber_lookup
        lookup = sample.create_group("fiber_lookup")
        lookup.create_dataset("fiber_ids", data=np.array([b"fiber1", b"fiber2"], dtype="S20"))
        lookup.create_dataset("chromosomes", data=np.array([b"chr1", b"chr1"], dtype="S10"))

        # chr1 group with required subgroups
        chrom = sample.create_group("chr1")
        chrom.create_dataset(
            "fiber_id_table", data=np.array([b"fiber1", b"fiber2"], dtype="S20")
        )

        fm = chrom.create_group("fiber_metadata")
        fm.create_dataset("fiber_int_ids", data=np.array([0, 1], dtype=np.uint32))
        fm.create_dataset("starts", data=np.array([100, 5000], dtype=np.uint32))
        fm.create_dataset("ends", data=np.array([2100, 8000], dtype=np.uint32))

    return db_path
