from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pytest

from PACKAGE.extract.qc import compute_qc_metrics, validate_bam


class _FakePopen:
    def __init__(self, lines):
        self.stdout = iter(lines)
        self.stderr = iter(())
        self.returncode = None

    def poll(self):
        return self.returncode

    def terminate(self):
        self.returncode = -15

    def wait(self, timeout=None):
        self.returncode = self.returncode if self.returncode is not None else 0
        return self.returncode


def _bam_and_index(tmp_path: Path) -> Path:
    bam = tmp_path / "sample.bam"
    bam.write_bytes(b"fake BAM")
    Path(f"{bam}.bai").write_bytes(b"fake index")
    return bam


def _mock_samtools(monkeypatch, *, sort_order="coordinate", tagged=True):
    def fake_run(command, **kwargs):
        if "quickcheck" in command:
            return SimpleNamespace(returncode=0, stdout="", stderr="")
        if "-H" in command:
            header = f"@HD\tVN:1.6\tSO:{sort_order}\n@SQ\tSN:chr1\tLN:1000\n"
            return SimpleNamespace(returncode=0, stdout=header, stderr="")
        raise AssertionError(command)

    tags = "\tMM:Z:C+m,0;\tML:B:C,200" if tagged else ""
    line = f"read1\t0\tchr1\t1\t60\t10M\t*\t0\t0\tACGT\tFFFF{tags}\n"
    monkeypatch.setattr("PACKAGE.extract.qc.subprocess.run", fake_run)
    monkeypatch.setattr(
        "PACKAGE.extract.qc.subprocess.Popen",
        lambda *args, **kwargs: _FakePopen([line]),
    )


def test_validate_bam_reports_header_and_tags(tmp_path, monkeypatch):
    bam = _bam_and_index(tmp_path)
    _mock_samtools(monkeypatch)
    result = validate_bam(bam, "ont", "/tools/samtools")
    assert result["sort_order"] == "coordinate"
    assert result["n_contigs"] == 1
    assert result["records_checked_for_tags"] == 1
    assert result["records_with_mm_ml_tags"] == 1


def test_validate_bam_requires_index(tmp_path):
    bam = tmp_path / "sample.bam"
    bam.write_bytes(b"fake BAM")
    with pytest.raises(FileNotFoundError, match="index not found"):
        validate_bam(bam, "ont")


def test_validate_bam_requires_coordinate_sort(tmp_path, monkeypatch):
    bam = _bam_and_index(tmp_path)
    _mock_samtools(monkeypatch, sort_order="queryname")
    with pytest.raises(ValueError, match="coordinate-sorted"):
        validate_bam(bam, "ont", "/tools/samtools")


def test_validate_ont_bam_requires_modification_tags(tmp_path, monkeypatch):
    bam = _bam_and_index(tmp_path)
    _mock_samtools(monkeypatch, tagged=False)
    with pytest.raises(ValueError, match="MM/ML"):
        validate_bam(bam, "ont", "/tools/samtools")


def test_compute_qc_metrics_lists_output_sizes(tmp_path):
    (tmp_path / "mods.tsv").write_text("abc")
    (tmp_path / "calls.bed").write_text("12345")
    result = compute_qc_metrics(tmp_path)
    assert result["files"]["mods.tsv"]["size_bytes"] == 3
    assert result["files"]["calls.bed"]["size_bytes"] == 5
