"""Tests for the database layer (read-only query API).

These tests use the small synthetic HDF5 from ``conftest.tiny_db_path``. As Batch 2
and Batch 3 of the port land, we add per-method tests here. Right now (post-Batch-1)
we test: open/close, metadata reading, sample/chromosome enumeration, fiber-count
summary, and the backward-compat annotation path resolution.
"""

from __future__ import annotations

from pathlib import Path

import h5py
import pandas as pd
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


# ===========================================================================
# Batch 2 tests: per-fiber data accessors
# ===========================================================================

# ---------------------------------------------------------------------------
# get_fibers_at — coordinate queries
# ---------------------------------------------------------------------------
def test_get_fibers_at_finds_overlapping_fibers(tiny_db_path_with_layers):
    """A region inside fiber1's span should return fiber1."""
    with FiberDatabase(tiny_db_path_with_layers) as db:
        fibers = db.get_fibers_at("chr1", 500, 1500, sample="d0")
    assert fibers == ["fiber1"]


def test_get_fibers_at_returns_multiple(tiny_db_path_with_layers):
    """A region spanning both fibers should return both."""
    with FiberDatabase(tiny_db_path_with_layers) as db:
        fibers = db.get_fibers_at("chr1", 0, 10_000, sample="d0")
    assert set(fibers) == {"fiber1", "fiber2"}


def test_get_fibers_at_returns_empty_for_no_overlap(tiny_db_path_with_layers):
    """A region with no fibers should return []."""
    with FiberDatabase(tiny_db_path_with_layers) as db:
        fibers = db.get_fibers_at("chr1", 3000, 4000, sample="d0")
    assert fibers == []


def test_get_fibers_at_strict_half_open_boundary(tiny_db_path_with_layers):
    """Region [end, end+N) should NOT include a fiber that ends at exactly `end`.

    fiber1 spans [100, 2100). A query for [2100, 3000) should return no fibers —
    fiber1's end coincides with the query start, no overlap. This is the half-open
    overlap semantics we agreed on (different from V8's inclusive comparison).
    """
    with FiberDatabase(tiny_db_path_with_layers) as db:
        fibers = db.get_fibers_at("chr1", 2100, 3000, sample="d0")
    assert fibers == []


def test_get_fibers_at_missing_chromosome(tiny_db_path_with_layers):
    """Querying a chromosome with no data should return []."""
    with FiberDatabase(tiny_db_path_with_layers) as db:
        fibers = db.get_fibers_at("chrZ", 0, 1000, sample="d0")
    assert fibers == []


# ---------------------------------------------------------------------------
# get_nucleosomes
# ---------------------------------------------------------------------------
def test_get_nucleosomes_returns_arrays(tiny_db_path_with_layers):
    """fiber1 has 3 nucleosomes in the fixture; check we get them all."""
    with FiberDatabase(tiny_db_path_with_layers) as db:
        nucs = db.get_nucleosomes("fiber1", "chr1", sample="d0")
    assert set(nucs.keys()) == {"starts", "ends", "widths", "linkers"}
    assert len(nucs["starts"]) == 3
    assert list(nucs["starts"]) == [150, 350, 550]
    assert list(nucs["widths"]) == [150, 150, 150]


def test_get_nucleosomes_unknown_fiber_returns_empty(tiny_db_path_with_layers):
    """Fiber not in the chromosome should return an empty dict, not error."""
    with FiberDatabase(tiny_db_path_with_layers) as db:
        nucs = db.get_nucleosomes("nonexistent_fiber", "chr1", sample="d0")
    assert nucs == {}


# ---------------------------------------------------------------------------
# get_methylation
# ---------------------------------------------------------------------------
def test_get_methylation_returns_all_fields(tiny_db_path_with_layers):
    """5mC fixture has positions, probabilities, is_methylated; we should get all."""
    with FiberDatabase(tiny_db_path_with_layers) as db:
        meth = db.get_methylation("fiber1", "chr1", mod_type="5mC", sample="d0")
    assert set(meth.keys()) == {"positions", "probabilities", "is_methylated"}
    assert len(meth["positions"]) == 4
    assert list(meth["positions"]) == [200, 400, 600, 800]
    assert list(meth["is_methylated"]) == [1, 0, 1, 0]


def test_get_methylation_missing_mod_type_returns_empty(tiny_db_path_with_layers):
    """6mA is not stored in this fixture; should get an empty dict."""
    with FiberDatabase(tiny_db_path_with_layers) as db:
        meth = db.get_methylation("fiber1", "chr1", mod_type="6mA", sample="d0")
    assert meth == {}


def test_get_methylation_rejects_invalid_mod_type(tiny_db_path_with_layers):
    """An unknown mod_type should raise ValueError from schema.layer_path."""
    # The slice-index path lookup will raise before we even check the data path.
    with FiberDatabase(tiny_db_path_with_layers) as db:
        with pytest.raises(ValueError, match="Unknown layer"):
            db.get_methylation("fiber1", "chr1", mod_type="not_a_real_mod", sample="d0")


# ---------------------------------------------------------------------------
# get_msp
# ---------------------------------------------------------------------------
def test_get_msp_returns_arrays(tiny_db_path_with_layers):
    """fiber1 has 2 MSPs in the fixture."""
    with FiberDatabase(tiny_db_path_with_layers) as db:
        msp = db.get_msp("fiber1", "chr1", sample="d0")
    assert set(msp.keys()) == {"starts", "ends", "widths"}
    assert len(msp["starts"]) == 2
    assert list(msp["starts"]) == [100, 600]


def test_get_msp_fiber_with_no_msps_returns_empty(tiny_db_path_with_layers):
    """A fiber not in the MSP slice index should get an empty dict."""
    # fiber1 and fiber2 both have MSPs in this fixture; query a fiber that doesn't exist.
    with FiberDatabase(tiny_db_path_with_layers) as db:
        msp = db.get_msp("nonexistent_fiber", "chr1", sample="d0")
    assert msp == {}


# ---------------------------------------------------------------------------
# _find_fiber_chromosome
# ---------------------------------------------------------------------------
def test_find_fiber_chromosome_returns_chrom(tiny_db_path_with_layers):
    """_find_fiber_chromosome looks up a fiber's chromosome from fiber_lookup."""
    with FiberDatabase(tiny_db_path_with_layers) as db:
        chrom = db._find_fiber_chromosome("fiber1", sample="d0")
    assert chrom == "chr1"


def test_find_fiber_chromosome_unknown_returns_none(tiny_db_path_with_layers):
    with FiberDatabase(tiny_db_path_with_layers) as db:
        assert db._find_fiber_chromosome("nonexistent_fiber", sample="d0") is None


# ===========================================================================
# Batch 3 tests: annotations + bulk queries + cross-sample
# ===========================================================================

# ---------------------------------------------------------------------------
# get_annotation_regions / get_annotation_by_id
# ---------------------------------------------------------------------------
def test_get_annotation_regions_returns_tuples(tiny_db_path_with_layers):
    """get_annotation_regions returns (chrom, start, end) tuples for the type."""
    with FiberDatabase(tiny_db_path_with_layers) as db:
        regions = db.get_annotation_regions("Typical_Enhancer")
    assert regions == [("chr1", 400, 700), ("chr1", 4900, 5500)]


def test_get_annotation_regions_unknown_returns_empty(tiny_db_path_with_layers):
    """Unknown annotation type should return []."""
    with FiberDatabase(tiny_db_path_with_layers) as db:
        regions = db.get_annotation_regions("Compartment_B")  # not in fixture
    assert regions == []


def test_get_annotation_by_id_returns_region(tiny_db_path_with_layers):
    """get_annotation_by_id strips trailing _N and finds by full ID."""
    with FiberDatabase(tiny_db_path_with_layers) as db:
        region = db.get_annotation_by_id("Typical_Enhancer_1")
    assert region == ("chr1", 400, 700)


def test_get_annotation_by_id_missing_returns_none(tiny_db_path_with_layers):
    with FiberDatabase(tiny_db_path_with_layers) as db:
        assert db.get_annotation_by_id("Typical_Enhancer_999") is None


# ---------------------------------------------------------------------------
# query_annotation_fast — the workhorse
# ---------------------------------------------------------------------------
def test_query_annotation_fast_returns_dataframe(tiny_db_path_with_layers):
    """Basic shape check: produces a DataFrame with expected columns."""
    with FiberDatabase(tiny_db_path_with_layers) as db:
        df = db.query_annotation_fast("Typical_Enhancer", sample="d0")
    assert isinstance(df, pd.DataFrame)
    assert not df.empty
    # Universal columns
    for col in ("region_id", "annotation", "fiber_id", "sample"):
        assert col in df.columns
    # Layer-specific columns (default is nucleosomes + 5mC)
    assert "n_nucleosomes" in df.columns
    assert "n_cpg" in df.columns


def test_query_annotation_fast_n_nucleosomes_correct(tiny_db_path_with_layers):
    """For Typical_Enhancer_1 region [400, 700), fiber1 has 2 nucleosomes fully inside.

    Fiber1 has nucleosomes at starts=[150, 350, 550] ends=[300, 500, 700].
    Within [400, 700) using V8's inclusive ``<=`` semantics:
      - n1 [150, 300]: start<400 → excluded
      - n2 [350, 500]: start>=400? No, 350<400 → excluded
      - n3 [550, 700]: start=550>=400 and end=700<=700 → included
    So fiber1 has 1 nucleosome in this region.
    """
    with FiberDatabase(tiny_db_path_with_layers) as db:
        df = db.query_annotation_fast(
            "Typical_Enhancer", sample="d0", feature_types=["nucleosomes"]
        )
    # Get the row for fiber1 at the first enhancer (region_id = "chr1:400-700")
    row = df[(df["fiber_id"] == "fiber1") & (df["region_id"] == "chr1:400-700")]
    assert len(row) == 1
    assert int(row["n_nucleosomes"].iloc[0]) == 1


def test_query_annotation_fast_pct_methylated_correct(tiny_db_path_with_layers):
    """In region [150, 850) (CGI), fiber1 has 4 CpGs: [200, 400, 600, 800].

    is_methylated values are [1, 0, 1, 0]. All four fall inside.
    Expected pct_methylated = (1+0+1+0)/4 * 100 = 50.0.
    """
    with FiberDatabase(tiny_db_path_with_layers) as db:
        df = db.query_annotation_fast(
            "CGI", sample="d0", feature_types=["5mC"]
        )
    row = df[df["fiber_id"] == "fiber1"]
    assert len(row) == 1
    assert int(row["n_cpg"].iloc[0]) == 4
    assert row["pct_methylated"].iloc[0] == 50.0


def test_query_annotation_fast_max_regions_caps_output(tiny_db_path_with_layers):
    """max_regions=1 should only process the first region."""
    with FiberDatabase(tiny_db_path_with_layers) as db:
        df_all = db.query_annotation_fast("Typical_Enhancer", sample="d0")
        df_capped = db.query_annotation_fast(
            "Typical_Enhancer", sample="d0", max_regions=1
        )
    # The capped version should have fewer rows (only Typical_Enhancer_1's results).
    assert df_capped["region_id"].nunique() == 1
    assert df_all["region_id"].nunique() == 2


def test_query_annotation_fast_unknown_annotation_returns_empty(
    tiny_db_path_with_layers,
):
    """Unknown annotation type yields an empty DataFrame, not an error."""
    with FiberDatabase(tiny_db_path_with_layers) as db:
        df = db.query_annotation_fast("Nonexistent_Type", sample="d0")
    assert df.empty


def test_query_annotation_alias_matches_fast(tiny_db_path_with_layers):
    """V8 had query_annotation as alias for query_annotation_fast; we preserve."""
    with FiberDatabase(tiny_db_path_with_layers) as db:
        df1 = db.query_annotation_fast("CGI", sample="d0")
        df2 = db.query_annotation("CGI", sample="d0")
    pd.testing.assert_frame_equal(df1, df2)


# ---------------------------------------------------------------------------
# query_by_id — single-region path
# ---------------------------------------------------------------------------
def test_query_by_id_returns_dataframe(tiny_db_path_with_layers):
    """query_by_id for a known annotation should produce a DataFrame."""
    with FiberDatabase(tiny_db_path_with_layers) as db:
        df = db.query_by_id("CGI_1", sample="d0")
    assert isinstance(df, pd.DataFrame)
    assert not df.empty
    assert df["region_id"].iloc[0] == "CGI_1"


def test_query_by_id_unknown_returns_empty(tiny_db_path_with_layers):
    with FiberDatabase(tiny_db_path_with_layers) as db:
        df = db.query_by_id("NonexistentID_999", sample="d0")
    assert df.empty


# ---------------------------------------------------------------------------
# compare_samples
# ---------------------------------------------------------------------------
def test_compare_samples_single_sample_matches_query(tiny_db_path_with_layers):
    """With one sample, compare_samples should equal query_annotation_fast."""
    with FiberDatabase(tiny_db_path_with_layers) as db:
        df_compare = db.compare_samples("CGI")
        df_query = db.query_annotation_fast("CGI", sample="d0")
    # compare_samples adds sample column at the end of pd.concat; reset_index for cmp.
    pd.testing.assert_frame_equal(
        df_compare.reset_index(drop=True), df_query.reset_index(drop=True)
    )


def test_compare_samples_empty_db_returns_empty(tmp_path):
    """Database with no annotations should produce empty result, not crash."""
    import h5py
    db_path = tmp_path / "empty_anno.h5"
    with h5py.File(db_path, "w") as f:
        f.create_group("d0")
    with FiberDatabase(db_path) as db:
        df = db.compare_samples("CGI")
    assert df.empty
