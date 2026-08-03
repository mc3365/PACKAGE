from __future__ import annotations

import csv
from pathlib import Path

import pandas as pd
import pytest

from PACKAGE.extract.pacbio import (
    convert_pacbio_all_5mc_to_modkit,
    inspect_pacbio_all,
    normalize_fibertools_block_bed,
)


_PACBIO_ALL = """\
#ct\tst\ten\tfiber\tscore\tstrand\tsam_flag\tHP\tRG\tfiber_length\tec\trq\ttotal_AT_bp\ttotal_m6a_bp\ttotal_nuc_bp\ttotal_msp_bp\ttotal_5mC_bp\tnuc_starts\tnuc_lengths\tref_nuc_starts\tref_nuc_lengths\tmsp_starts\tmsp_lengths\tfire\tref_msp_starts\tref_msp_lengths\tm6a\tref_m6a\tm6a_qual\t5mC\tref_5mC\t5mC_qual
chr1\t100\t500\tfiber_a\t8\t+\t0\tUNK\trg1\t400\t0\t1\t0\t2\t2\t2\t3\t0,100\t50,80\t100,200\t50,80\t0,120\t30,40\t0,244\t100,220\t30,40\t5,10\t105,110\t255,125\t7,9,12\t107,-1,112\t255,200,128
chr1\t700\t900\tfiber_b\t4\t-\t0\tUNK\trg1\t200\t0\t1\t0\t1\t1\t1\t1\t0\t80\t700\t80\t0\t60\t0\t720\t60\t3\t703\t250\t5\t705\t64
"""


def test_inspect_pacbio_all_validates_named_columns_and_content(tmp_path: Path):
    path = tmp_path / "pacbio_all.tsv"
    path.write_text(_PACBIO_ALL)

    summary = inspect_pacbio_all(path)

    assert summary["records_checked"] == 2
    assert summary["has_fire_columns"] is True
    assert summary["total_5mC_calls_checked"] == 4
    assert summary["mapped_5mC_calls_checked"] == 3
    assert summary["positive_fire_entries_checked"] == 1


def test_inspect_pacbio_all_rejects_missing_semantic_columns(tmp_path: Path):
    path = tmp_path / "bad_all.tsv"
    path.write_text("#ct\tst\ten\tfiber\tstrand\tref_5mC\nchr1\t0\t1\tread\t+\t10\n")

    with pytest.raises(ValueError, match="missing required PacBio 5mC columns"):
        inspect_pacbio_all(path)


def test_convert_pacbio_all_5mc_to_modkit_compatible_tsv(tmp_path: Path):
    all_path = tmp_path / "pacbio_all.tsv"
    out_path = tmp_path / "pacbio_5mc_modkit.tsv"
    all_path.write_text(_PACBIO_ALL)

    summary = convert_pacbio_all_5mc_to_modkit(all_path, out_path)

    assert summary["converted_5mC_rows"] == 3
    assert summary["skipped_unmapped_5mC_rows"] == 1
    with out_path.open() as f:
        rows = list(csv.reader(f, delimiter="\t"))
    assert rows[0][0] == "read_id"
    assert rows[1][0] == "fiber_a"
    assert rows[1][2] == "107"
    assert rows[1][3] == "chr1"
    assert rows[1][12] == "1.000000"
    assert rows[1][13] == "m"
    assert rows[-1][0] == "fiber_b"
    assert rows[-1][12] == f"{64 / 255:.6f}"


def test_pacbio_all_missing_dot_positions_are_skipped(tmp_path: Path):
    all_path = tmp_path / "pacbio_all_dot.tsv"
    out_path = tmp_path / "pacbio_5mc_modkit.tsv"
    all_path.write_text(
        _PACBIO_ALL
        + "chr2\t1000\t1200\tfiber_c\t3\t+\t0\tUNK\trg1\t200\t0\t1\t0\t1\t0\t0\t1\t"
        + ".\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t42\t.\t.\n"
    )

    summary = inspect_pacbio_all(all_path)
    converted = convert_pacbio_all_5mc_to_modkit(all_path, out_path)

    assert summary["total_5mC_calls_checked"] == 5
    assert summary["mapped_5mC_calls_checked"] == 3
    assert converted["converted_5mC_rows"] == 3
    assert converted["skipped_unmapped_5mC_rows"] == 2


def test_normalize_fibertools_block_bed_accepts_bed12(tmp_path: Path):
    raw = tmp_path / "msp.raw.bed"
    out = tmp_path / "msp.bed"
    raw.write_text(
        "chr1\t100\t500\tfiber_a\t8\t+\t100\t500\t147,112,219\t"
        "4\t0,50,80,1\t0,40,220,399\n"
    )

    summary = normalize_fibertools_block_bed(raw, out, expected="interval")

    assert summary["rows"] == 1
    assert out.read_text().strip().split("\t")[9] == "4"


def test_normalize_fibertools_block_bed_converts_bed9_plus_2(tmp_path: Path):
    raw = tmp_path / "6ma.raw.bed"
    out = tmp_path / "6ma.bed"
    raw.write_text(
        "chr1\t100\t500\tfiber_a\t8\t+\t100\t500\t128,0,12247\t"
        "0,1,1,1\t0,10,20,399\n"
    )

    summary = normalize_fibertools_block_bed(raw, out, expected="single_base")
    cols = out.read_text().strip().split("\t")

    assert summary["single_base_blocks"] == 2
    assert len(cols) == 12
    assert cols[8] == "128,0,12247"
    assert cols[9] == "4"
    assert cols[10] == "0,1,1,1"
    assert cols[11] == "0,10,20,399"


def test_normalize_fibertools_block_bed_converts_bed8_plus_2(tmp_path: Path):
    raw = tmp_path / "5mc.raw.bed"
    out = tmp_path / "5mc.bed"
    raw.write_text(
        "chr1\t100\t500\tfiber_a\t8\t+\t100\t500\t"
        "139,69,10,1\t0,100,220,399\n"
    )

    normalize_fibertools_block_bed(raw, out)
    cols = out.read_text().strip().split("\t")

    assert len(cols) == 12
    assert cols[8] == "0"
    assert cols[9] == "4"
    assert cols[10] == "139,69,10,1"
    assert cols[11] == "0,100,220,399"


def test_normalize_fibertools_block_bed_rejects_bad_block_lists(tmp_path: Path):
    raw = tmp_path / "bad.bed"
    raw.write_text(
        "chr1\t100\t500\tfiber_a\t8\t+\t100\t500\t128,0,12247\t"
        "0,1,1\t0,10\n"
    )

    with pytest.raises(ValueError, match="Could not identify consistent block columns"):
        normalize_fibertools_block_bed(raw, tmp_path / "bad.out.bed")


def test_converted_pacbio_5mc_can_be_read_like_modkit(tmp_path: Path):
    all_path = tmp_path / "pacbio_all.tsv"
    out_path = tmp_path / "pacbio_5mc_modkit.tsv"
    all_path.write_text(_PACBIO_ALL)
    convert_pacbio_all_5mc_to_modkit(all_path, out_path)

    df = pd.read_csv(out_path, sep="\t")

    assert df["mod_code"].tolist() == ["m", "m", "m"]
    assert df["ref_position"].tolist() == [107, 112, 705]
    assert df["mod_qual"].between(0, 1).all()
