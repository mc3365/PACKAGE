"""Tests for PACKAGE.db.builder.

Batch 1 scope: foundation + format parsers. Builds tiny in-memory text fixtures
representing real ft-extract / modkit-extract output and verifies each parser
produces the right structure.

Subsequent batches will add:
  - Per-chromosome assembly tests (Batch 2)
  - End-to-end build test against a synthetic database (Batch 3)
"""

from __future__ import annotations

from pathlib import Path

import pytest

from PACKAGE.config import Config
from PACKAGE.db.builder import FiberDatabaseBuilder, build_database

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

    monkeypatch.setattr("PACKAGE.db.builder.FiberDatabaseBuilder", MockBuilder)
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
