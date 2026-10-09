"""Tests for the schema module.

Since the schema is the contract between builder and reader, these tests verify the
invariants that both sides depend on. If any of these fail, we have a schema bug.
"""

from __future__ import annotations

import pytest

from mei_fiber.db import schema


def test_supported_layers_have_field_definitions():
    """Every layer in SUPPORTED_LAYERS must have an entry in LAYER_FIELDS."""
    for layer in schema.SUPPORTED_LAYERS:
        assert layer in schema.LAYER_FIELDS, (
            f"Layer {layer!r} listed in SUPPORTED_LAYERS but missing from LAYER_FIELDS"
        )


def test_layer_path_rejects_unknown():
    """layer_path() should raise on layers not in SUPPORTED_LAYERS."""
    with pytest.raises(ValueError, match="Unknown layer"):
        schema.layer_path("d0", "chr1", "not_a_real_layer")


def test_slice_index_path_rejects_unknown():
    """slice_index_path() should raise on layers not in SUPPORTED_LAYERS."""
    with pytest.raises(ValueError, match="Unknown layer"):
        schema.slice_index_path("d0", "chr1", "not_a_real_layer")


def test_path_helpers_produce_consistent_structure():
    """Path helpers should produce paths consistent with the layout in the docstring."""
    # All sample-level paths should start with the sample name
    assert schema.fiber_lookup_path("mysample").startswith("mysample/")
    assert schema.chrom_group_path("mysample", "chr1") == "mysample/chr1"
    assert schema.layer_path("mysample", "chr1", "5mC") == "mysample/chr1/5mC"
    # The slice-index path should be under _indices/
    assert "_indices" in schema.slice_index_path("mysample", "chr1", "5mC")


def test_reserved_top_level_not_a_sample():
    """RESERVED_TOP_LEVEL names should not be confused with sample names anywhere."""
    # Schema only declares them; this just guards against accidentally adding a
    # sample-only path helper that conflicts.
    assert "annotations" in schema.RESERVED_TOP_LEVEL
    assert "metadata" in schema.RESERVED_TOP_LEVEL


def test_slice_index_dtype_fields():
    """SLICE_INDEX_DTYPE must have the three fields readers/writers expect."""
    assert schema.SLICE_INDEX_DTYPE.names == ("fiber_int_id", "start", "end")


def test_annotation_region_dtype_fields():
    """ANNOTATION_REGION_DTYPE must have the four expected fields."""
    assert schema.ANNOTATION_REGION_DTYPE.names == ("chr", "start", "end", "id")


def test_schema_version_is_string():
    """SCHEMA_VERSION must be a string (gets written as HDF5 attr)."""
    assert isinstance(schema.SCHEMA_VERSION, str)
    assert "." in schema.SCHEMA_VERSION  # at least major.minor
