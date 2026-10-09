"""Tests for mei_fiber.db.builder.

Batch 1 scope: foundation + format parsers. Builds tiny in-memory text fixtures
representing real ft-extract / modkit-extract output and verifies each parser
produces the right structure.

Subsequent batches will add:
  - Per-chromosome assembly tests (Batch 2)
  - End-to-end build test against a synthetic database (Batch 3)
"""

from __future__ import annotations

from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import pytest

from mei_fiber.config import Config
from mei_fiber.db.builder import FiberDatabaseBuilder, build_database
from mei_fiber.db.database import FiberDatabase

# ---------------------------------------------------------------------------
# Tiny fixture inputs
# ---------------------------------------------------------------------------
# These mimic the formats produced by modkit / fibertools-rs ft extract. Each is
# kept small enough to hand-verify what the parser should produce.

# BED12 with single-base blocks (6mA format).
# Cols: chrom, start, end, name(read_id), score, strand, thickStart, thickEnd,
#       rgb, blockCount, blockSizes (csv), blockStarts (csv)
# First and last blocks are boundary markers (the parser skips them).
_BED12_6MA_TEXT = """\
chr1\t1000\t2000\tfiber_a\t0\t+\t1000\t2000\t128,0,128\t5\t1,1,1,1,1\t0,100,200,300,1000
chr1\t1000\t2000\tfiber_b\t0\t+\t1000\t2000\t128,0,128\t4\t1,1,1,1\t0,50,150,1000
chr2\t5000\t6000\tfiber_c\t0\t+\t5000\t6000\t128,0,128\t3\t1,1,1\t0,500,1000
"""

# BED12 with interval blocks (MSP format). Same column layout, sizes can vary.
# We want the parser to filter blocks smaller than min_size, and skip first/last.
_BED12_MSP_TEXT = """\
chr1\t1000\t2000\tfiber_a\t0\t+\t1000\t2000\t147,112,219\t4\t1,50,8,1\t0,100,500,1000
chr1\t1000\t2000\tfiber_b\t0\t+\t1000\t2000\t147,112,219\t3\t1,80,1\t0,200,1000
"""

_FIRE_ACCESSIBILITY_TEXT = """\
chr1\t110\t170\tfiber_a\t100\t+\t110\t170\t147,112,219\t0.05\tH1
chr1\t500\t760\tfiber_b\t101\t+\t500\t760\t169,169,169\t1.01\tUNK
chr1\t900\t950\tfiber_unknown\t101\t+\t900\t950\t169,169,169\t1.01\tUNK
chr2\t100\t150\tfiber_c\t101\t+\t100\t150\t169,169,169\t0.25\tH2
"""

# modkit `extract full` TSV — first line is header (we skip it).
# Real header has ~20 columns but our parser only uses 5 of them (cols 0,2,3,12,13).
# We pad with dummies to reach >= 14 cols.
def _modkit_row(read_id: str, ref_pos: int, chrom: str, prob: float, code: str) -> str:
    cols = [""] * 14
    cols[0] = read_id
    cols[2] = str(ref_pos)
    cols[3] = chrom
    cols[12] = str(prob)
    cols[13] = code
    return "\t".join(cols)


_MODKIT_TEXT = "\n".join([
    "read_id\tfwd_pos\tref_pos\tchrom\t...\t...\t...\t...\t...\t...\t...\t...\tmod_qual\tmod_code",  # header
    _modkit_row("fiber_a", 100, "chr1", 0.95, "m"),  # 5mC, kept
    _modkit_row("fiber_a", 100, "chr1", 0.10, "h"),  # 5hmC, kept
    _modkit_row("fiber_a", 200, "chr1", 0.30, "m"),  # 5mC, low prob — still kept
    _modkit_row("fiber_b", 50, "chr1", 0.85, "m"),
    _modkit_row("fiber_unknown", 999, "chr1", 0.5, "m"),  # NOT in known_fibers -> dropped
    _modkit_row("fiber_a", 0, "chrZ", 0.5, "m"),         # chrZ not in known_fibers -> dropped
]) + "\n"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _write(tmp_path: Path, name: str, text: str) -> Path:
    p = tmp_path / name
    p.write_text(text)
    return p


def _make_builder(tmp_path: Path) -> FiberDatabaseBuilder:
    """Construct a builder with a minimal config (no real files needed for parser tests)."""
    cfg = Config.model_validate({
        "output_dir": str(tmp_path / "out"),
        "samples": [
            {
                "name": "test",
                "layers": {"nucleosomes": str(tmp_path / "n.csv")},
            },
        ],
    })
    return FiberDatabaseBuilder(cfg, [cfg.samples[0]])


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------
def test_build_database_calls_builder(tmp_path, monkeypatch):
    """build_database should resolve samples and dispatch to the builder."""
    # Use monkeypatch so we don't need a real implementation of build() to test
    # the dispatch logic.
    called = {}

    class MockBuilder(FiberDatabaseBuilder):
        def build(self):
            called["yes"] = True

    monkeypatch.setattr("mei_fiber.db.builder.FiberDatabaseBuilder", MockBuilder)
    cfg = Config.model_validate({
        "output_dir": str(tmp_path),
        "samples": [
            {"name": "d0", "layers": {"nucleosomes": str(tmp_path / "n0.csv")}},
            {"name": "d4", "layers": {"nucleosomes": str(tmp_path / "n4.csv")}},
        ],
    })
    build_database(cfg, samples=["d0"])
    assert called.get("yes") is True


def test_build_database_unknown_sample_raises(tmp_path):
    """Asking to build a sample not in the config should raise."""
    cfg = Config.model_validate({
        "output_dir": str(tmp_path),
        "samples": [{"name": "d0", "layers": {"nucleosomes": str(tmp_path / "n.csv")}}],
    })
    with pytest.raises(KeyError):
        build_database(cfg, samples=["notreal"])


# ---------------------------------------------------------------------------
# validate_inputs
# ---------------------------------------------------------------------------
def test_validate_inputs_missing_file_lists_what_is_missing(tmp_path):
    """validate_inputs should report every missing file in one error."""
    cfg = Config.model_validate({
        "output_dir": str(tmp_path),
        "samples": [
            {
                "name": "d0",
                "layers": {
                    "nucleosomes": str(tmp_path / "missing_nuc.csv"),
                    "5mC":         str(tmp_path / "missing_mods.tsv.gz"),
                },
            },
        ],
    })
    builder = FiberDatabaseBuilder(cfg, list(cfg.samples))
    with pytest.raises(FileNotFoundError) as exc:
        builder.validate_inputs()
    assert "missing_nuc.csv" in str(exc.value)
    assert "missing_mods.tsv.gz" in str(exc.value)


def test_validate_inputs_all_present_no_error(tmp_path):
    """If all configured files exist, validate_inputs returns silently."""
    n_path = _write(tmp_path, "n.csv", "chrom,read_id,nuc_start,nuc_end,width,linker_dist\n")
    cfg = Config.model_validate({
        "output_dir": str(tmp_path),
        "samples": [{"name": "d0", "layers": {"nucleosomes": str(n_path)}}],
    })
    builder = FiberDatabaseBuilder(cfg, list(cfg.samples))
    builder.validate_inputs()  # must not raise


# ---------------------------------------------------------------------------
# parse_bed12_single_base (6mA format)
# ---------------------------------------------------------------------------
def test_parse_bed12_single_base_skips_boundaries_and_unknowns(tmp_path):
    builder = _make_builder(tmp_path)
    bed_path = _write(tmp_path, "6ma.bed", _BED12_6MA_TEXT)
    # known_fibers includes fiber_a and fiber_b on chr1, but only fiber_c on chr2.
    known = {
        "chr1": {b"fiber_a", b"fiber_b"},
        "chr2": {b"fiber_c"},
    }
    result = builder.parse_bed12_single_base(bed_path, known)

    # fiber_a has 5 blocks; we skip first + last, leaving 3 positions:
    # abs_start=1000, offsets[1:-1] = [100, 200, 300] -> positions = [1100, 1200, 1300]
    assert result["chr1"]["fids"].count(b"fiber_a") == 3
    fiber_a_positions = [
        p for p, f in zip(result["chr1"]["pos"], result["chr1"]["fids"], strict=True)
        if f == b"fiber_a"
    ]
    assert fiber_a_positions == [1100, 1200, 1300]

    # fiber_b has 4 blocks; skip first + last, leaves 2: offsets [50, 150] -> [1050, 1150]
    fiber_b_positions = [
        p for p, f in zip(result["chr1"]["pos"], result["chr1"]["fids"], strict=True)
        if f == b"fiber_b"
    ]
    assert fiber_b_positions == [1050, 1150]

    # fiber_c on chr2 has 3 blocks; skip first + last leaves only 1 (offset=500)
    assert result["chr2"]["fids"] == [b"fiber_c"]
    assert result["chr2"]["pos"] == [5500]


def test_parse_bed12_single_base_drops_unknown_fiber(tmp_path):
    """Reads not in known_fibers must be dropped."""
    builder = _make_builder(tmp_path)
    bed_path = _write(tmp_path, "6ma.bed", _BED12_6MA_TEXT)
    # Empty known set for chr1
    known: dict[str, set[bytes]] = {"chr1": set(), "chr2": {b"fiber_c"}}
    result = builder.parse_bed12_single_base(bed_path, known)
    # chr1 should have no entries; chr2 should still have fiber_c
    assert "chr1" not in result or len(result["chr1"]["fids"]) == 0
    assert b"fiber_c" in result["chr2"]["fids"]


def test_parse_bed12_single_base_drops_unknown_chrom(tmp_path):
    """Rows on chroms not in known_fibers must be dropped."""
    builder = _make_builder(tmp_path)
    bed_path = _write(tmp_path, "6ma.bed", _BED12_6MA_TEXT)
    # known_fibers only has chr1
    known = {"chr1": {b"fiber_a", b"fiber_b"}}
    result = builder.parse_bed12_single_base(bed_path, known)
    # chr2 must not appear in output even though there's a row in input
    assert "chr2" not in result


# ---------------------------------------------------------------------------
# parse_bed12_intervals (MSP format)
# ---------------------------------------------------------------------------
def test_parse_bed12_intervals_skips_small_blocks_and_boundaries(tmp_path):
    builder = _make_builder(tmp_path)
    bed_path = _write(tmp_path, "msp.bed", _BED12_MSP_TEXT)
    known = {"chr1": {b"fiber_a", b"fiber_b"}}
    result = builder.parse_bed12_intervals(bed_path, known, min_size=10)

    # fiber_a: 4 blocks [1, 50, 8, 1]. Skip first+last -> [50, 8].
    # min_size=10 filters out the 8 -> only the 50-bp block survives.
    fiber_a_widths = [
        w for w, f in zip(result["chr1"]["widths"], result["chr1"]["fids"], strict=True)
        if f == b"fiber_a"
    ]
    assert fiber_a_widths == [50]
    # That block's start = abs_start(1000) + offset(100) = 1100; end = 1100 + 50 = 1150
    fiber_a_starts = [
        s for s, f in zip(result["chr1"]["starts"], result["chr1"]["fids"], strict=True)
        if f == b"fiber_a"
    ]
    fiber_a_ends = [
        e for e, f in zip(result["chr1"]["ends"], result["chr1"]["fids"], strict=True)
        if f == b"fiber_a"
    ]
    assert fiber_a_starts == [1100]
    assert fiber_a_ends == [1150]

    # fiber_b: 3 blocks [1, 80, 1]. Skip first+last -> [80]. 80 >= 10 -> kept.
    fiber_b_widths = [
        w for w, f in zip(result["chr1"]["widths"], result["chr1"]["fids"], strict=True)
        if f == b"fiber_b"
    ]
    assert fiber_b_widths == [80]


def test_parse_bed12_intervals_no_min_size_keeps_all(tmp_path):
    """min_size=0 should keep every non-boundary block."""
    builder = _make_builder(tmp_path)
    bed_path = _write(tmp_path, "msp.bed", _BED12_MSP_TEXT)
    known = {"chr1": {b"fiber_a", b"fiber_b"}}
    result = builder.parse_bed12_intervals(bed_path, known, min_size=0)
    # fiber_a's full middle two blocks: 50 and 8 -> both kept
    fiber_a_widths = sorted(
        w for w, f in zip(result["chr1"]["widths"], result["chr1"]["fids"], strict=True)
        if f == b"fiber_a"
    )
    assert fiber_a_widths == [8, 50]


# ---------------------------------------------------------------------------
# parse_fire_accessibility_bed (ft fire --extract format)
# ---------------------------------------------------------------------------
def test_parse_fire_accessibility_bed_uses_column_10_scores(tmp_path):
    builder = _make_builder(tmp_path)
    fire_path = _write(tmp_path, "acc.model.results.sort.bed", _FIRE_ACCESSIBILITY_TEXT)
    known = {"chr1": {b"fiber_a", b"fiber_b"}}

    result = builder.parse_fire_accessibility_bed(fire_path, known)

    assert result["chr1"]["fids"] == [b"fiber_a", b"fiber_b"]
    assert result["chr1"]["starts"] == [110, 500]
    assert result["chr1"]["ends"] == [170, 760]
    assert result["chr1"]["widths"] == [60, 260]
    assert result["chr1"]["scores"] == [0.05, 1.01]
    assert result["chr1"]["haplotypes"] == [b"H1", b"UNK"]
    assert "chr2" not in result


# ---------------------------------------------------------------------------
# parse_modkit_extract (5mC + 5hmC)
# ---------------------------------------------------------------------------
def test_parse_modkit_extract_routes_to_both_layers(tmp_path):
    builder = _make_builder(tmp_path)
    mods_path = _write(tmp_path, "raw_mods.tsv", _MODKIT_TEXT)
    known = {"chr1": {b"fiber_a", b"fiber_b"}}
    mc, hmc = builder.parse_modkit_extract(mods_path, known)

    # Expected from fixture:
    #   fiber_a (chr1): 5mC at 100 (prob 0.95), 5mC at 200 (prob 0.30), 5hmC at 100 (0.10)
    #   fiber_b (chr1): 5mC at 50 (prob 0.85)
    #   fiber_unknown: dropped (not in known_fibers)
    #   chrZ row: dropped (chrom not in known_fibers)
    assert len(mc["chr1"]["pos"]) == 3
    assert sorted(mc["chr1"]["pos"]) == [50, 100, 200]
    assert len(hmc["chr1"]["pos"]) == 1
    assert hmc["chr1"]["pos"] == [100]


def test_parse_modkit_extract_skips_layer_if_not_wanted(tmp_path):
    """want_5hmc=False should leave hmc empty."""
    builder = _make_builder(tmp_path)
    mods_path = _write(tmp_path, "raw_mods.tsv", _MODKIT_TEXT)
    known = {"chr1": {b"fiber_a", b"fiber_b"}}
    mc, hmc = builder.parse_modkit_extract(mods_path, known, want_5hmc=False)
    # 5mC should still be populated
    assert len(mc["chr1"]["pos"]) == 3
    # 5hmC should be empty
    assert len(hmc) == 0


def test_parse_modkit_extract_unknown_chrom_dropped(tmp_path):
    """Rows with chrom not in known_fibers must be dropped."""
    builder = _make_builder(tmp_path)
    mods_path = _write(tmp_path, "raw_mods.tsv", _MODKIT_TEXT)
    known = {"chr1": {b"fiber_a"}}  # only fiber_a on chr1
    mc, _ = builder.parse_modkit_extract(mods_path, known)
    # fiber_b row on chr1 should be dropped (not in known set)
    # fiber_a's two 5mC calls should remain
    assert b"fiber_b" not in mc.get("chr1", {}).get("fids", [])
    assert mc["chr1"]["fids"].count(b"fiber_a") == 2


# ---------------------------------------------------------------------------
# Class-level constants (sanity check)
# ---------------------------------------------------------------------------
def test_builder_class_constants_match_v8():
    """Column indices for modkit parsing should match V8's layout."""
    assert FiberDatabaseBuilder._MODKIT_COL_READ_ID == 0
    assert FiberDatabaseBuilder._MODKIT_COL_REF_POS == 2
    assert FiberDatabaseBuilder._MODKIT_COL_CHROM == 3
    assert FiberDatabaseBuilder._MODKIT_COL_MOD_QUAL == 12
    assert FiberDatabaseBuilder._MODKIT_COL_MOD_CODE == 13


# ---------------------------------------------------------------------------
# Batch 2: integer ID mapping + per-chromosome assembly/write helpers
# ---------------------------------------------------------------------------
def test_build_fiber_int_map_normalizes_to_bytes_and_uint32(tmp_path):
    builder = _make_builder(tmp_path)
    str_to_int, table = builder._build_fiber_int_map(["fiber_b", b"fiber_a"])

    assert table.dtype.kind == "S"
    assert table.tolist() == [b"fiber_b", b"fiber_a"]
    assert str_to_int[b"fiber_b"] == np.uint32(0)
    assert str_to_int[b"fiber_a"] == np.uint32(1)

    mapped = builder._map_fids_to_int(["fiber_a", b"fiber_b"], str_to_int)
    assert mapped.dtype == np.uint32
    assert mapped.tolist() == [1, 0]


def test_sort_and_assemble_stable_sorts_parallel_arrays(tmp_path):
    builder = _make_builder(tmp_path)
    str_to_int, _ = builder._build_fiber_int_map(["fiber_a", "fiber_b"])

    sorted_fids, sorted_arrays = builder._sort_and_assemble(
        ["fiber_b", "fiber_a", "fiber_b", "fiber_a"],
        {"pos": [40, 10, 50, 20], "prob": [0.4, 0.1, 0.5, 0.2]},
        str_to_int,
    )

    assert sorted_fids.tolist() == [0, 0, 1, 1]
    # Stable sort preserves within-fiber input order: fiber_a positions 10 then 20.
    assert sorted_arrays["pos"].tolist() == [10, 20, 40, 50]
    assert sorted_arrays["prob"].tolist() == [0.1, 0.2, 0.4, 0.5]


def test_load_nucleosome_inputs_builds_fiber_metadata_and_known_sets(tmp_path):
    nuc_path = _write(
        tmp_path,
        "nuc.csv",
        "\n".join([
            "chrom,read_id,nuc_start,nuc_end,width,linker_dist",
            "chr1,fiber_b,300,450,150,20",
            "chr1,fiber_a,100,250,150,",
            "chr1,fiber_a,260,400,140,10",
            "chr2,fiber_c,500,650,150,.",
        ]).replace(",.", ",") + "\n",
    )
    cfg = Config.model_validate({
        "output_dir": str(tmp_path / "out"),
        "samples": [{"name": "d0", "layers": {"nucleosomes": str(nuc_path)}}],
    })
    builder = FiberDatabaseBuilder(cfg, list(cfg.samples))

    nuc_df, fiber_meta, known = builder.load_nucleosome_inputs(cfg.samples[0])

    assert len(nuc_df) == 4
    assert set(known["chr1"]) == {b"fiber_a", b"fiber_b"}
    row = fiber_meta[fiber_meta["read_id"] == "fiber_a"].iloc[0]
    assert row["start"] == 100
    assert row["end"] == 400


def _assembly_inputs():
    nuc_df = pd.DataFrame({
        "chrom": ["chr1", "chr1", "chr1", "chr1", "chr1"],
        "read_id": ["fiber_b", "fiber_a", "fiber_a", "fiber_b", "fiber_a"],
        "nuc_start": [500, 100, 300, 700, 200],
        "nuc_end": [650, 180, 360, 820, 280],
        "width": [150, 80, 60, 120, 80],
        "linker_dist": [-1, -1, 20, 50, 20],
    })
    fiber_meta = (
        nuc_df.groupby(["chrom", "read_id"])
        .agg(start=("nuc_start", "min"), end=("nuc_end", "max"))
        .reset_index()
        .sort_values(["chrom", "read_id"])
        .reset_index(drop=True)
    )
    mc_data = {
        "chr1": {
            "fids": [b"fiber_b", b"fiber_a", b"fiber_b"],
            "pos": [710, 120, 730],
            "prob": [0.49, 0.50, 0.90],
        }
    }
    hmc_data = {"chr1": {"fids": [b"fiber_a"], "pos": [150], "prob": [0.20]}}
    ma_data = {"chr1": {"fids": [b"fiber_b", b"fiber_a"], "pos": [740, 130]}}
    msp_data = {"chr1": {"fids": [b"fiber_b"], "starts": [705], "ends": [760], "widths": [55]}}
    fire_data = {
        "chr1": {
            "fids": [b"fiber_b", b"fiber_a"],
            "starts": [700, 110],
            "ends": [760, 170],
            "widths": [60, 60],
            "scores": [1.01, 0.05],
            "haplotypes": [b"UNK", b"H1"],
        }
    }
    return nuc_df, fiber_meta, mc_data, hmc_data, ma_data, msp_data, fire_data


def test_assemble_chromosome_data_matches_v8_sorting_and_threshold(tmp_path):
    builder = _make_builder(tmp_path)
    nuc_df, fiber_meta, mc_data, hmc_data, ma_data, msp_data, fire_data = _assembly_inputs()

    data = builder.assemble_chromosome_data(
        "chr1", fiber_meta, nuc_df, mc_data, hmc_data, ma_data, msp_data, fire_data
    )

    # Sorted unique read IDs define per-chromosome integer IDs.
    assert data["fiber_id_table"].tolist() == [b"fiber_a", b"fiber_b"]
    assert data["fm_int_ids"].tolist() == [0, 1]

    # Nucleosomes are stable-sorted by int ID. fiber_a rows appear in original
    # fiber_a order from nuc_df: starts 100, 300, 200.
    assert data["nuc_int_ids"].tolist() == [0, 0, 0, 1, 1]
    assert data["nuc_starts"].tolist() == [100, 300, 200, 500, 700]

    # Methylation calls are grouped by fiber int ID; threshold is >= 0.5.
    assert data["mc_int_ids"].tolist() == [0, 1, 1]
    assert data["mc_pos"].tolist() == [120, 710, 730]
    assert data["mc_is_meth"].tolist() == [True, False, True]

    assert data["hmc_int_ids"].tolist() == [0]
    assert data["ma_int_ids"].tolist() == [0, 1]
    assert data["msp_int_ids"].tolist() == [1]
    assert data["fire_int_ids"].tolist() == [0, 1]
    assert data["fire_scores"].tolist() == pytest.approx([0.05, 1.01])
    assert data["fire_haplotypes"].tolist() == [b"H1", b"UNK"]


def test_write_chrom_data_and_indices_round_trip_through_fiber_database(tmp_path):
    builder = _make_builder(tmp_path)
    nuc_df, fiber_meta, mc_data, hmc_data, ma_data, msp_data, fire_data = _assembly_inputs()
    data = builder.assemble_chromosome_data(
        "chr1", fiber_meta, nuc_df, mc_data, hmc_data, ma_data, msp_data, fire_data
    )

    db_path = tmp_path / "assembled.h5"
    with h5py.File(db_path, "w") as hf:
        meta = hf.create_group("metadata")
        meta.attrs["version"] = "8.0"
        meta.attrs["genome_version"] = "mm10"
        meta.attrs["methylation_threshold"] = 0.5

        sample = hf.create_group("d0")
        builder._write_chrom_data(sample, "chr1", data)
        builder._build_indices(sample, "chr1", data)

        lookup = sample.create_group("fiber_lookup")
        lookup.create_dataset(
            "fiber_ids",
            data=np.array([b"fiber_a", b"fiber_b"], dtype="S50"),
            compression="gzip",
        )
        lookup.create_dataset(
            "chromosomes",
            data=np.array([b"chr1", b"chr1"], dtype="S10"),
            compression="gzip",
        )

    with h5py.File(db_path, "r") as hf:
        assert hf["d0/chr1/_indices/nucleosomes_slices"][:].tolist() == [
            (0, 0, 3),
            (1, 3, 5),
        ]
        assert hf["d0/chr1/_indices/5mC_slices"][:].tolist() == [
            (0, 0, 1),
            (1, 1, 3),
        ]
        assert "5hmC" in hf["d0/chr1"]
        assert "6mA" in hf["d0/chr1"]
        assert "msp" in hf["d0/chr1"]
        assert "fire_accessibility" in hf["d0/chr1"]
        assert hf["d0/chr1/_indices/fire_accessibility_slices"][:].tolist() == [
            (0, 0, 1),
            (1, 1, 2),
        ]

    with FiberDatabase(db_path) as db:
        nuc = db.get_nucleosomes("fiber_a", "chr1", sample="d0")
        assert nuc["starts"].tolist() == [100, 300, 200]

        meth = db.get_methylation("fiber_b", "chr1", sample="d0")
        assert meth["positions"].tolist() == [710, 730]
        assert meth["is_methylated"].tolist() == [False, True]

        hmc = db.get_methylation("fiber_a", "chr1", mod_type="5hmC", sample="d0")
        assert hmc["positions"].tolist() == [150]

        msp = db.get_msp("fiber_b", "chr1", sample="d0")
        assert msp["starts"].tolist() == [705]

        fire = db.get_fire_accessibility("fiber_a", "chr1", sample="d0")
        assert fire["starts"].tolist() == [110]
        assert fire["scores"].tolist() == pytest.approx([0.05])


def test_write_chrom_data_skips_absent_optional_layers(tmp_path):
    builder = _make_builder(tmp_path)
    nuc_df, fiber_meta, *_ = _assembly_inputs()
    data = builder.assemble_chromosome_data("chr1", fiber_meta, nuc_df)

    db_path = tmp_path / "nuc_only.h5"
    with h5py.File(db_path, "w") as hf:
        sample = hf.create_group("d0")
        builder._write_chrom_data(sample, "chr1", data)
        builder._build_indices(sample, "chr1", data)

        assert "nucleosomes" in sample["chr1"]
        assert "nucleosomes_slices" in sample["chr1/_indices"]
        assert "5mC" not in sample["chr1"]
        assert "5hmC" not in sample["chr1"]
        assert "6mA" not in sample["chr1"]
        assert "msp" not in sample["chr1"]


def test_build_database_end_to_end_from_extracted_files(tmp_path):
    """Tiny end-to-end build from extracted-style inputs, then query the result."""
    nuc_path = _write(
        tmp_path,
        "nuc.csv",
        "\n".join([
            "chrom,read_id,nuc_start,nuc_end,width,linker_dist",
            "chr1,fiber_a,100,180,80,",
            "chr1,fiber_a,220,300,80,40",
            "chr1,fiber_b,500,650,150,",
            "chr2,fiber_c,1000,1150,150,",
        ]) + "\n",
    )
    mods_path = _write(tmp_path, "raw_mods.tsv", _MODKIT_TEXT)
    ma_path = _write(tmp_path, "6ma.bed", _BED12_6MA_TEXT)
    msp_path = _write(tmp_path, "msp.bed", _BED12_MSP_TEXT)
    fire_path = _write(tmp_path, "fire.bed", _FIRE_ACCESSIBILITY_TEXT)
    anno_path = _write(
        tmp_path,
        "master.bed",
        "\n".join([
            "chr1\t90\t320\tCGI_1",
            "chr1\t450\t700\tCGI_2",
            "chr2\t950\t1200\tEnhancer_1",
        ]) + "\n",
    )

    cfg = Config.model_validate({
        "output_dir": str(tmp_path / "out"),
        "output_file": "built.h5",
        "samples": [
            {
                "name": "d0",
                "layers": {
                    "nucleosomes": str(nuc_path),
                    "5mC": str(mods_path),
                    "5hmC": str(mods_path),
                    "6mA": str(ma_path),
                    "msp": str(msp_path),
                    "fire_accessibility": str(fire_path),
                },
            }
        ],
        "annotations": {"master": str(anno_path)},
        "parameters": {"methylation_threshold": 0.5, "min_msp_size": 10},
    })

    build_database(cfg)

    assert cfg.output_path.exists()
    with h5py.File(cfg.output_path, "r") as hf:
        assert "d0/chr1/fiber_id_table" in hf
        assert "d0/chr1/_indices/nucleosomes_slices" in hf
        assert "d0/chr1/_indices/fire_accessibility_slices" in hf
        assert "d0/fiber_lookup/fiber_ids" in hf
        assert "annotations/master/features/CGI" in hf
        assert hf["metadata"].attrs["methylation_threshold"] == 0.5
        assert hf["metadata"].attrs["mei_fiber_version"] == "0.5.0"
        assert hf["metadata"].attrs["package_version"] == "0.5.0"

    with FiberDatabase(cfg.output_path) as db:
        assert db.samples == ["d0"]
        assert db.list_annotations() == ["CGI", "Enhancer"]
        assert db.get_chromosomes("d0") == ["chr1", "chr2"]

        fiber_a_nuc = db.get_nucleosomes("fiber_a", "chr1", sample="d0")
        assert fiber_a_nuc["starts"].tolist() == [100, 220]

        fiber_a_mc = db.get_methylation("fiber_a", "chr1", sample="d0")
        assert fiber_a_mc["positions"].tolist() == [100, 200]
        assert fiber_a_mc["is_methylated"].tolist() == [True, False]

        df = db.query_annotation_fast("CGI", sample="d0", max_regions=1)
        assert set(df["fiber_id"]) == {"fiber_a"}
        assert df.iloc[0]["n_nucleosomes"] == 2

        fire_df = db.query_annotation_fast(
            "CGI",
            sample="d0",
            feature_types=["fire_accessibility"],
            max_regions=1,
        )
        assert fire_df.iloc[0]["n_fire_accessibility"] == 1
        assert fire_df.iloc[0]["min_fire_score"] == pytest.approx(0.05)
