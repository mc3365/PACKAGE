"""Shared pytest fixtures.

Fixtures defined here are auto-discovered by pytest — no import needed in test files.

This file provides synthetic HDF5 databases for testing the query layer without
needing real Fiber-seq data on disk. As batches land, the fixtures grow:

  - ``tiny_db_path``: minimal (Batch 1) — metadata + fiber_metadata only
  - ``tiny_db_path_with_layers``: tiny_db_path PLUS layer data (Batch 2) —
    nucleosomes, 5mC, msp populated for 2 fibers on chr1
"""

from __future__ import annotations

from pathlib import Path

import h5py
import numpy as np
import pytest


@pytest.fixture
def tiny_db_path(tmp_path: Path) -> Path:
    """Minimal PACKAGE-shaped HDF5 in a per-test tmpdir.

    Contains top-level metadata and one (sample, chromosome) with fiber_metadata
    for 2 fibers. No modification layers. Used by Batch 1 tests.
    """
    db_path = tmp_path / "tiny.h5"
    with h5py.File(db_path, "w") as f:
        meta = f.create_group("metadata")
        meta.attrs["version"] = "8.0"
        meta.attrs["genome_version"] = "mm10"
        meta.attrs["methylation_threshold"] = 0.5

        sample = f.create_group("d0")
        lookup = sample.create_group("fiber_lookup")
        lookup.create_dataset(
            "fiber_ids", data=np.array([b"fiber1", b"fiber2"], dtype="S20")
        )
        lookup.create_dataset(
            "chromosomes", data=np.array([b"chr1", b"chr1"], dtype="S10")
        )

        chrom = sample.create_group("chr1")
        chrom.create_dataset(
            "fiber_id_table", data=np.array([b"fiber1", b"fiber2"], dtype="S20")
        )
        fm = chrom.create_group("fiber_metadata")
        fm.create_dataset("fiber_int_ids", data=np.array([0, 1], dtype=np.uint32))
        fm.create_dataset("starts", data=np.array([100, 5000], dtype=np.uint32))
        fm.create_dataset("ends", data=np.array([2100, 8000], dtype=np.uint32))

    return db_path


# ---------------------------------------------------------------------------
# Layer-populated fixture for Batch 2
# ---------------------------------------------------------------------------

# Define the data layout once so the fixture is easy to read. Each fiber gets:
#  - a set of nucleosome intervals (starts, ends, widths, linkers)
#  - a set of 5mC calls (positions, probabilities, is_methylated)
#  - a set of MSP intervals (starts, ends, widths)
#
# Slice indices encode the row range per fiber within each layer's data arrays.

_NUC_DATA = {
    # fiber0: 3 nucleosomes at rows [0:3]; fiber1: 2 nucleosomes at rows [3:5]
    "fiber_int_ids": np.array([0, 0, 0, 1, 1], dtype=np.uint32),
    "starts":        np.array([150, 350, 550, 5100, 5300], dtype=np.uint32),
    "ends":          np.array([300, 500, 700, 5250, 5450], dtype=np.uint32),
    "widths":        np.array([150, 150, 150, 150, 150], dtype=np.uint32),
    "linkers":       np.array([-1, 50, 50, -1, 50], dtype=np.int32),
}
_NUC_SLICES = np.array(
    [(0, 0, 3), (1, 3, 5)],
    dtype=[("fiber_int_id", np.uint32), ("start", np.uint32), ("end", np.uint32)],
)

_M5C_DATA = {
    # fiber0: 4 CpGs at rows [0:4]; fiber1: 2 CpGs at rows [4:6]
    "fiber_int_ids":  np.array([0, 0, 0, 0, 1, 1], dtype=np.uint32),
    "positions":      np.array([200, 400, 600, 800, 5200, 5400], dtype=np.uint32),
    "probabilities":  np.array([0.9, 0.1, 0.8, 0.2, 0.95, 0.3], dtype=np.float32),
    "is_methylated":  np.array([1, 0, 1, 0, 1, 0], dtype=np.uint8),
}
_M5C_SLICES = np.array(
    [(0, 0, 4), (1, 4, 6)],
    dtype=[("fiber_int_id", np.uint32), ("start", np.uint32), ("end", np.uint32)],
)

_MSP_DATA = {
    # fiber0: 2 MSPs at rows [0:2]; fiber1: 1 MSP at rows [2:3]
    "fiber_int_ids": np.array([0, 0, 1], dtype=np.uint32),
    "starts":        np.array([100, 600, 5100], dtype=np.uint32),
    "ends":          np.array([200, 750, 5400], dtype=np.uint32),
    "widths":        np.array([100, 150, 300], dtype=np.uint32),
}
_MSP_SLICES = np.array(
    [(0, 0, 2), (1, 2, 3)],
    dtype=[("fiber_int_id", np.uint32), ("start", np.uint32), ("end", np.uint32)],
)

_FIRE_DATA = {
    # fiber0: two FIRE/accessibility calls; fiber1: one call
    "fiber_int_ids": np.array([0, 0, 1], dtype=np.uint32),
    "starts": np.array([120, 650, 5200], dtype=np.uint32),
    "ends": np.array([180, 760, 5300], dtype=np.uint32),
    "widths": np.array([60, 110, 100], dtype=np.uint32),
    "scores": np.array([0.05, 1.01, 0.25], dtype=np.float64),
    "haplotypes": np.array([b"H1", b"H1", b"UNK"], dtype="S10"),
}
_FIRE_SLICES = np.array(
    [(0, 0, 2), (1, 2, 3)],
    dtype=[("fiber_int_id", np.uint32), ("start", np.uint32), ("end", np.uint32)],
)


@pytest.fixture
def tiny_db_path_with_layers(tmp_path: Path) -> Path:
    """PACKAGE-shaped HDF5 with all the structure Batch 2 methods need.

    Two fibers on chr1:
      - fiber1 (int_id=0): spans [100, 2100), has 3 nucleosomes, 4 5mC calls, 2 MSPs
      - fiber2 (int_id=1): spans [5000, 8000), has 2 nucleosomes, 2 5mC calls, 1 MSP

    Positions chosen so a region query like [400, 700) selects only fiber1's middle
    features, exercising the per-fiber slice logic.
    """
    db_path = tmp_path / "tiny_layers.h5"
    with h5py.File(db_path, "w") as f:
        # Top-level metadata
        meta = f.create_group("metadata")
        meta.attrs["version"] = "8.0"
        meta.attrs["genome_version"] = "mm10"
        meta.attrs["methylation_threshold"] = 0.5

        # Sample + fiber_lookup
        sample = f.create_group("d0")
        lookup = sample.create_group("fiber_lookup")
        lookup.create_dataset(
            "fiber_ids", data=np.array([b"fiber1", b"fiber2"], dtype="S20")
        )
        lookup.create_dataset(
            "chromosomes", data=np.array([b"chr1", b"chr1"], dtype="S10")
        )

        # chr1 group
        chrom = sample.create_group("chr1")
        chrom.create_dataset(
            "fiber_id_table", data=np.array([b"fiber1", b"fiber2"], dtype="S20")
        )
        fm = chrom.create_group("fiber_metadata")
        fm.create_dataset("fiber_int_ids", data=np.array([0, 1], dtype=np.uint32))
        fm.create_dataset("starts", data=np.array([100, 5000], dtype=np.uint32))
        fm.create_dataset("ends", data=np.array([2100, 8000], dtype=np.uint32))

        # Layer data + slice indices
        indices = chrom.create_group("_indices")
        for layer_name, layer_data, slices in [
            ("nucleosomes", _NUC_DATA, _NUC_SLICES),
            ("5mC",         _M5C_DATA, _M5C_SLICES),
            ("msp",         _MSP_DATA, _MSP_SLICES),
            ("fire_accessibility", _FIRE_DATA, _FIRE_SLICES),
        ]:
            layer_grp = chrom.create_group(layer_name)
            for field, arr in layer_data.items():
                layer_grp.create_dataset(field, data=arr)
            indices.create_dataset(f"{layer_name}_slices", data=slices)

        # Annotation features: two annotation types.
        # Region [400, 700) overlaps fiber1's middle features but not fiber2.
        # Region [4900, 5500) overlaps fiber2's MSPs and some methylation.
        anno_dtype = [("chr", "S10"), ("start", np.uint32),
                      ("end", np.uint32), ("id", "S50")]
        enhancer_data = np.array([
            (b"chr1", 400, 700, b"Typical_Enhancer_1"),
            (b"chr1", 4900, 5500, b"Typical_Enhancer_2"),
        ], dtype=anno_dtype)
        cgi_data = np.array([
            (b"chr1", 150, 850, b"CGI_1"),
        ], dtype=anno_dtype)
        anno_features = f.create_group("annotations/master/features")
        anno_features.create_dataset("Typical_Enhancer", data=enhancer_data)
        anno_features.create_dataset("CGI", data=cgi_data)

    return db_path


@pytest.fixture
def tiny_db_with_monster_fiber(tmp_path: Path) -> Path:
    """A db with 2 normal fibers + 1 monster fiber (span > 100 kb).

    Used to verify FiberDatabase.MAX_FIBER_SPAN filtering across all overlap
    code paths (spatial index, array scan, bulk query).

    Layout:
      - fiber_normal_a: int_id=0, spans [100, 2100)        — normal, 2 kb
      - fiber_monster:  int_id=1, spans [500, 500_000_000)  — monster, way too large
      - fiber_normal_b: int_id=2, spans [5000, 8000)       — normal, 3 kb

    A region query for [1000, 2000) should match only fiber_normal_a, NOT the
    monster (even though the monster's span trivially overlaps everything).
    """
    db_path = tmp_path / "tiny_with_monster.h5"
    with h5py.File(db_path, "w") as f:
        meta = f.create_group("metadata")
        meta.attrs["version"] = "8.0"
        meta.attrs["genome_version"] = "mm10"

        sample = f.create_group("d0")
        lookup = sample.create_group("fiber_lookup")
        lookup.create_dataset(
            "fiber_ids",
            data=np.array([b"fiber_normal_a", b"fiber_monster", b"fiber_normal_b"], dtype="S30"),
        )
        lookup.create_dataset(
            "chromosomes", data=np.array([b"chr1", b"chr1", b"chr1"], dtype="S10")
        )

        chrom = sample.create_group("chr1")
        chrom.create_dataset(
            "fiber_id_table",
            data=np.array([b"fiber_normal_a", b"fiber_monster", b"fiber_normal_b"], dtype="S30"),
        )
        fm = chrom.create_group("fiber_metadata")
        fm.create_dataset("fiber_int_ids", data=np.array([0, 1, 2], dtype=np.uint32))
        fm.create_dataset("starts", data=np.array([100, 500, 5000], dtype=np.uint32))
        fm.create_dataset("ends",   data=np.array([2100, 500_000_000, 8000], dtype=np.uint32))

    return db_path
