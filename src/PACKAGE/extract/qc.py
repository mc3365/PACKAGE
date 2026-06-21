"""Cheap pre-extraction validation and output inventory metrics."""

from __future__ import annotations

import subprocess
from pathlib import Path
from typing import Any


def _bam_index_path(bam_path: Path) -> Path | None:
    candidates = [Path(f"{bam_path}.bai"), bam_path.with_suffix(".bai")]
    return next((path for path in candidates if path.exists()), None)


def _sample_modification_tags(
    bam_path: Path,
    samtools_executable: str,
    max_records: int,
) -> tuple[int, int]:
    """Return (records_checked, records_with_paired_MM_ML_tags)."""
    process = subprocess.Popen(
        [samtools_executable, "view", str(bam_path)],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    checked = 0
    tagged = 0
    assert process.stdout is not None
    try:
        for line in process.stdout:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 11:
                continue
            checked += 1
            tags = {field[:2] for field in fields[11:] if len(field) >= 5}
            if ({"MM", "ML"} <= tags) or ({"Mm", "Ml"} <= tags):
                tagged += 1
            if checked >= max_records:
                break
    finally:
        if process.poll() is None:
            process.terminate()
        process.wait(timeout=10)
    return checked, tagged


def validate_bam(
    bam_path: Path,
    platform: str,
    samtools_executable: str = "samtools",
    max_records: int = 1_000,
) -> dict[str, Any]:
    """Validate readability, index, sort order, contigs, and modification tags."""
    if not bam_path.exists():
        raise FileNotFoundError(f"BAM not found: {bam_path}")
    if not bam_path.is_file():
        raise ValueError(f"BAM path is not a file: {bam_path}")
    index_path = _bam_index_path(bam_path)
    if index_path is None:
        raise FileNotFoundError(
            f"BAM index not found for {bam_path}; expected {bam_path}.bai or "
            f"{bam_path.with_suffix('.bai')}"
        )

    quickcheck = subprocess.run(
        [samtools_executable, "quickcheck", "-v", str(bam_path)],
        check=False,
        capture_output=True,
        text=True,
    )
    if quickcheck.returncode != 0:
        detail = quickcheck.stderr.strip() or quickcheck.stdout.strip()
        raise ValueError(f"samtools quickcheck failed for {bam_path}: {detail}")

    header_result = subprocess.run(
        [samtools_executable, "view", "-H", str(bam_path)],
        check=True,
        capture_output=True,
        text=True,
    )
    header_lines = header_result.stdout.splitlines()
    hd = next((line for line in header_lines if line.startswith("@HD")), "")
    sort_order = next(
        (field[3:] for field in hd.split("\t") if field.startswith("SO:")),
        "unknown",
    )
    if sort_order != "coordinate":
        raise ValueError(
            f"BAM must be coordinate-sorted; header reports SO:{sort_order}"
        )
    contigs = [
        field[3:]
        for line in header_lines
        if line.startswith("@SQ")
        for field in line.split("\t")
        if field.startswith("SN:")
    ]
    if not contigs:
        raise ValueError("BAM header contains no @SQ reference records")
    if not any(contig.startswith("chr") for contig in contigs):
        raise ValueError(
            "BAM reference names do not use the expected 'chr' prefix convention"
        )

    platform = platform.lower()
    checked, tagged = _sample_modification_tags(
        bam_path, samtools_executable, max_records
    )
    if checked == 0:
        raise ValueError(f"BAM contains no aligned records: {bam_path}")
    if platform == "ont" and tagged == 0:
        raise ValueError(
            f"No paired MM/ML modification tags found in the first {checked} records"
        )

    return {
        "bam": str(bam_path),
        "bam_size_bytes": bam_path.stat().st_size,
        "index": str(index_path),
        "sort_order": sort_order,
        "n_contigs": len(contigs),
        "records_checked_for_tags": checked,
        "records_with_mm_ml_tags": tagged,
    }


def compute_qc_metrics(sample_dir: Path) -> dict[str, Any]:
    """Return a cheap inventory of files produced in a sample output directory."""
    files = {}
    if sample_dir.exists():
        for path in sorted(sample_dir.iterdir()):
            if path.is_file():
                files[path.name] = {"path": str(path), "size_bytes": path.stat().st_size}
    return {"output_directory": str(sample_dir), "files": files}
