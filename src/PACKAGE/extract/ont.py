"""ONT extraction: aligned, modification-tagged BAM to builder inputs.

This module orchestrates the same tools used for the validated d0/d4 database:

* ``modkit extract full --cpg`` for the shared 5mC/5hmC table
* ``ft extract --m6a`` for 6mA BED12
* ``ft extract --msp`` for MSP BED12
* ``ft extract --nuc`` followed by streaming BED12-to-CSV conversion

Dorado basecalling and alignment remain upstream prerequisites. The input BAM must
already contain the modification tags consumed by modkit and fibertools-rs.
"""

from __future__ import annotations

import csv
import shutil
import subprocess
from collections.abc import Sequence
from pathlib import Path

from PACKAGE.config import Config, SampleConfig
from PACKAGE.extract.manifest import capture_tool_versions, write_manifest
from PACKAGE.extract.qc import compute_qc_metrics, validate_bam
from PACKAGE.utils import get_logger, smart_open

log = get_logger(__name__)


def _resolve_executable(name: str) -> str:
    """Return an executable path or raise a clear dependency error."""
    resolved = shutil.which(name)
    if resolved is None:
        raise FileNotFoundError(
            f"Required executable {name!r} was not found on PATH. "
            "Activate the ONT environment or set its executable in config.extraction."
        )
    return resolved


def _nucleosome_bed_path(csv_path: Path) -> Path:
    """Intermediate BED12 path paired with a nucleosome feature CSV."""
    return csv_path.with_suffix(".bed12")


def flatten_nucleosome_bed12(input_bed: Path, output_csv: Path) -> int:
    """Stream ft nucleosome BED12 into the CSV expected by the database builder.

    Linker distance is measured from the previous nucleosome end on the same BED12
    record. Zero-width boundary anchors are ignored, matching the validated legacy
    converter without retaining hundreds of millions of rows in memory.

    Returns:
        Number of nucleosome rows written.
    """
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    rows_written = 0
    with smart_open(input_bed) as source, output_csv.open("w", newline="") as target:
        writer = csv.writer(target)
        writer.writerow(
            ["chrom", "read_id", "nuc_start", "nuc_end", "width", "linker_dist"]
        )
        for line_number, line in enumerate(source, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 12:
                raise ValueError(
                    f"Invalid BED12 row at {input_bed}:{line_number}: "
                    f"expected at least 12 columns, found {len(cols)}"
                )
            chrom = cols[0]
            abs_start = int(cols[1])
            read_id = cols[3]
            sizes = [int(x) for x in cols[10].rstrip(",").split(",") if x]
            offsets = [int(x) for x in cols[11].rstrip(",").split(",") if x]
            if len(sizes) != len(offsets):
                raise ValueError(
                    f"Invalid BED12 blocks at {input_bed}:{line_number}: "
                    f"{len(sizes)} sizes != {len(offsets)} offsets"
                )

            previous_end: int | None = None
            for size, offset in zip(sizes, offsets, strict=True):
                if size == 0:
                    continue
                start = abs_start + offset
                end = start + size
                linker = "" if previous_end is None else start - previous_end
                writer.writerow([chrom, read_id, start, end, size, linker])
                previous_end = end
                rows_written += 1
    return rows_written


def _run(command: Sequence[str], command_log: list[list[str]] | None = None) -> None:
    """Run one external command with logging and fail on nonzero exit."""
    normalized = [str(part) for part in command]
    log.info("Running: " + " ".join(normalized))
    if command_log is not None:
        command_log.append(normalized)
    subprocess.run(normalized, check=True)


def _validate_outputs(sample: SampleConfig) -> dict[str, dict[str, str | int]]:
    """Require every configured extraction output to exist and be nonempty."""
    outputs: dict[str, dict[str, str | int]] = {}
    for layer, path in sample.layers.items():
        if not path.exists():
            raise FileNotFoundError(f"Configured {layer} output was not created: {path}")
        size = path.stat().st_size
        if size == 0:
            raise ValueError(f"Configured {layer} output is empty: {path}")
        outputs[layer] = {"path": str(path), "size_bytes": size}
    return outputs


def _validate_sample_inputs(config: Config, sample: SampleConfig) -> None:
    if sample.bam is None:
        raise ValueError(f"Sample {sample.name!r} has no BAM path in config")
    if not sample.bam.exists():
        raise FileNotFoundError(f"BAM not found for sample {sample.name!r}: {sample.bam}")
    if config.reference is None:
        raise ValueError("config.reference is required for ONT modkit extraction")
    if not config.reference.exists():
        raise FileNotFoundError(f"Reference FASTA not found: {config.reference}")

    methylation_paths = {
        path for layer, path in sample.layers.items() if layer in {"5mC", "5hmC"}
    }
    if len(methylation_paths) > 1:
        raise ValueError(
            f"Sample {sample.name!r} config points 5mC and 5hmC to different files. "
            "modkit extract produces both channels in one shared TSV; use one path."
        )


def extract_sample(config: Config, sample: SampleConfig) -> None:
    """Extract all configured ONT layers for one sample."""
    _validate_sample_inputs(config, sample)
    settings = config.extraction
    ft = _resolve_executable(settings.ft_executable)
    assert sample.bam is not None
    assert config.reference is not None
    command_log: list[list[str]] = []

    bam_qc = {}
    samtools: str | None = None
    if settings.validate_bam:
        samtools = _resolve_executable(settings.samtools_executable)
        log.info(f"Validating BAM: {sample.bam}")
        bam_qc = validate_bam(sample.bam, "ont", samtools)
        log.info(
            "BAM validation passed: "
            f"{bam_qc['records_with_mm_ml_tags']}/"
            f"{bam_qc['records_checked_for_tags']} sampled reads have MM/ML tags"
        )

    for path in sample.layers.values():
        path.parent.mkdir(parents=True, exist_ok=True)

    methylation_paths = {
        path for layer, path in sample.layers.items() if layer in {"5mC", "5hmC"}
    }
    if methylation_paths:
        modkit = _resolve_executable(settings.modkit_executable)
        output = next(iter(methylation_paths))
        if settings.overwrite or not output.exists():
            _run(
                [
                    modkit,
                    "extract",
                    "full",
                    "--cpg",
                    "--threads",
                    str(settings.threads),
                    "--reference",
                    str(config.reference),
                    str(sample.bam),
                    str(output),
                ],
                command_log,
            )
        else:
            log.info(f"Skipping existing modkit output: {output}")

    ft_layers = (("6mA", "--m6a"), ("msp", "--msp"))
    for layer, flag in ft_layers:
        output = sample.layers.get(layer)
        if output is None:
            continue
        if settings.overwrite or not output.exists():
            _run([ft, "extract", flag, str(output), str(sample.bam)], command_log)
        else:
            log.info(f"Skipping existing {layer} output: {output}")

    nuc_csv = sample.layers["nucleosomes"]
    nuc_bed = _nucleosome_bed_path(nuc_csv)
    if settings.overwrite or not nuc_csv.exists():
        if settings.overwrite or not nuc_bed.exists():
            _run(
                [ft, "extract", "--nuc", str(nuc_bed), str(sample.bam)],
                command_log,
            )
        rows = flatten_nucleosome_bed12(nuc_bed, nuc_csv)
        log.info(f"Wrote {rows:,} nucleosome rows: {nuc_csv}")
        if not settings.keep_nucleosome_bed:
            nuc_bed.unlink()
    else:
        log.info(f"Skipping existing nucleosome CSV: {nuc_csv}")

    outputs = _validate_outputs(sample)
    if settings.write_manifest:
        executables = {"fibertools-rs": ft}
        if methylation_paths:
            executables["modkit"] = modkit
        if samtools is not None:
            executables["samtools"] = samtools
        output_dir = nuc_csv.parent
        manifest_path = write_manifest(
            output_dir,
            {
                "sample": sample.name,
                "platform": "ont",
                "input": {
                    "bam": str(sample.bam),
                    "bam_size_bytes": sample.bam.stat().st_size,
                    "reference": str(config.reference),
                },
                "bam_qc": bam_qc,
                "tools": capture_tool_versions(executables),
                "parameters": {
                    "threads": settings.threads,
                    "overwrite": settings.overwrite,
                    "methylation_threshold": config.parameters.methylation_threshold,
                },
                "commands_run": command_log,
                "outputs": outputs,
                "output_inventory": compute_qc_metrics(output_dir),
            },
            filename=f"PACKAGE_manifest_{sample.name}.json",
        )
        log.info(f"Wrote extraction manifest: {manifest_path}")


def extract_samples(config: Config, samples: list[str] | None = None) -> None:
    """Run the ONT extraction pipeline for selected configured samples."""
    sample_names = samples or config.sample_names
    log.info(f"ONT extraction: samples={sample_names}")
    for name in sample_names:
        extract_sample(config, config.get_sample(name))
