from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from PACKAGE.config import Config
from PACKAGE.extract.ont import extract_sample, flatten_nucleosome_bed12

_NUC_BED12 = """\
chr1\t100\t500\tfiber_a\t0\t+\t100\t500\t0\t4\t0,80,70,0\t0,20,150,400
chr2\t1000\t1400\tfiber_b\t0\t-\t1000\t1400\t0\t3\t0,100,0\t0,50,400
"""


def _config(tmp_path: Path, *, overwrite: bool = False) -> Config:
    bam = tmp_path / "sample.bam"
    ref = tmp_path / "genome.fa"
    bam.write_bytes(b"BAM fixture")
    ref.write_text(">chr1\nACGT\n")
    return Config.model_validate(
        {
            "output_dir": str(tmp_path / "db"),
            "reference": str(ref),
            "samples": [
                {
                    "name": "d0",
                    "bam": str(bam),
                    "layers": {
                        "nucleosomes": str(tmp_path / "out" / "d0_nuc_features.csv"),
                        "5mC": str(tmp_path / "out" / "d0_raw_mods.tsv.gz"),
                        "5hmC": str(tmp_path / "out" / "d0_raw_mods.tsv.gz"),
                        "6mA": str(tmp_path / "out" / "d0_6ma.bed"),
                        "msp": str(tmp_path / "out" / "d0_msp.bed"),
                    },
                }
            ],
            "extraction": {"threads": 7, "overwrite": overwrite},
        }
    )


def test_flatten_nucleosome_bed12_streams_expected_features(tmp_path):
    bed = tmp_path / "nucleosomes.bed12"
    csv_path = tmp_path / "nucleosomes.csv"
    bed.write_text(_NUC_BED12)

    count = flatten_nucleosome_bed12(bed, csv_path)
    df = pd.read_csv(csv_path)

    assert count == 3
    assert df.columns.tolist() == [
        "chrom",
        "read_id",
        "nuc_start",
        "nuc_end",
        "width",
        "linker_dist",
    ]
    assert df["nuc_start"].tolist() == [120, 250, 1050]
    assert df["nuc_end"].tolist() == [200, 320, 1150]
    assert pd.isna(df.iloc[0]["linker_dist"])
    assert df.iloc[1]["linker_dist"] == 50
    assert pd.isna(df.iloc[2]["linker_dist"])


def test_flatten_nucleosome_bed12_rejects_mismatched_blocks(tmp_path):
    bed = tmp_path / "bad.bed12"
    bed.write_text("chr1\t0\t10\tread\t0\t+\t0\t10\t0\t2\t1,2\t0\n")
    with pytest.raises(ValueError, match="sizes !="):
        flatten_nucleosome_bed12(bed, tmp_path / "bad.csv")


def test_extract_sample_runs_validated_ont_commands(tmp_path, monkeypatch):
    cfg = _config(tmp_path)
    sample = cfg.samples[0]
    commands: list[list[str]] = []

    monkeypatch.setattr(
        "PACKAGE.extract.ont.shutil.which",
        lambda name: f"/tools/{name}",
    )

    def fake_run(command, check):
        assert check is True
        commands.append(command)
        if "--nuc" in command:
            Path(command[3]).write_text(_NUC_BED12)

    monkeypatch.setattr("PACKAGE.extract.ont.subprocess.run", fake_run)

    extract_sample(cfg, sample)

    assert commands[0] == [
        "/tools/modkit",
        "extract",
        "full",
        "--cpg",
        "--threads",
        "7",
        "--reference",
        str(cfg.reference),
        str(sample.bam),
        str(sample.layers["5mC"]),
    ]
    assert commands[1] == [
        "/tools/ft",
        "extract",
        "--m6a",
        str(sample.layers["6mA"]),
        str(sample.bam),
    ]
    assert commands[2] == [
        "/tools/ft",
        "extract",
        "--msp",
        str(sample.layers["msp"]),
        str(sample.bam),
    ]
    assert commands[3][0:3] == ["/tools/ft", "extract", "--nuc"]
    assert sample.layers["nucleosomes"].exists()


def test_extract_sample_skips_existing_outputs(tmp_path, monkeypatch):
    cfg = _config(tmp_path)
    sample = cfg.samples[0]
    for path in set(sample.layers.values()):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("existing")

    monkeypatch.setattr(
        "PACKAGE.extract.ont.shutil.which",
        lambda name: f"/tools/{name}",
    )
    commands = []
    monkeypatch.setattr(
        "PACKAGE.extract.ont.subprocess.run",
        lambda command, check: commands.append(command),
    )

    extract_sample(cfg, sample)
    assert commands == []


def test_extract_sample_requires_bam_and_reference(tmp_path):
    cfg = _config(tmp_path)
    cfg.samples[0].bam.unlink()
    with pytest.raises(FileNotFoundError, match="BAM not found"):
        extract_sample(cfg, cfg.samples[0])


def test_nucleosome_only_extraction_does_not_require_modkit(tmp_path, monkeypatch):
    cfg = _config(tmp_path)
    sample = cfg.samples[0]
    sample.layers = {"nucleosomes": sample.layers["nucleosomes"]}

    def fake_which(name):
        return "/tools/ft" if name == "ft" else None

    monkeypatch.setattr("PACKAGE.extract.ont.shutil.which", fake_which)

    def fake_run(command, check):
        assert check is True
        Path(command[3]).write_text(_NUC_BED12)

    monkeypatch.setattr("PACKAGE.extract.ont.subprocess.run", fake_run)
    extract_sample(cfg, sample)
    assert sample.layers["nucleosomes"].exists()
