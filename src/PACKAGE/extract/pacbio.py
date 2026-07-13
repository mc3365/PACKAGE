"""PacBio extraction and normalization.

PacBio Fiber-seq BAMs are handled differently from ONT at the extraction stage,
but the database builder should still receive the same normalized layer files:

* nucleosomes: PACKAGE nucleosome CSV
* 5mC: modkit-compatible TSV with probabilities in column 13
* 6mA: BED12-style single-base blocks
* MSP: BED12-style interval blocks

The important guardrail here is content-aware validation.  PacBio fibertools
outputs can be BED12, BED9+2, or BED8+2 depending on the track.  Rather than
trusting a fixed column number, this module inspects the fields by meaning:
coordinate columns, block size/start lists, named ``ft extract --all`` columns,
and per-row list lengths.
"""

from __future__ import annotations

import csv
import shutil
import subprocess
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

from PACKAGE.config import Config, SampleConfig
from PACKAGE.extract.manifest import capture_tool_versions, write_manifest
from PACKAGE.extract.ont import flatten_nucleosome_bed12
from PACKAGE.extract.qc import compute_qc_metrics, validate_bam
from PACKAGE.utils import get_logger, smart_open

log = get_logger(__name__)


_ALL_REQUIRED_5MC_COLUMNS = {"ct", "st", "en", "fiber", "strand", "ref_5mC", "5mC_qual"}
_ALL_FIRE_COLUMNS = {"ref_msp_starts", "ref_msp_lengths", "fire"}


@dataclass(frozen=True)
class BlockBedRecord:
    """Parsed fibertools BED-like row with normalized block lists."""

    columns: list[str]
    block_count_index: int | None
    sizes_index: int
    starts_index: int
    sizes: list[int]
    starts: list[int]


def _resolve_executable(name: str) -> str:
    resolved = shutil.which(name)
    if resolved is None:
        raise FileNotFoundError(
            f"Required executable {name!r} was not found on PATH. "
            "Activate the PacBio/Fiber-seq environment or set it in config.extraction."
        )
    return resolved


def _run(command: Sequence[str], command_log: list[list[str]] | None = None) -> None:
    normalized = [str(part) for part in command]
    log.info("Running: " + " ".join(normalized))
    if command_log is not None:
        command_log.append(normalized)
    subprocess.run(normalized, check=True)


def _insert_before_suffix(path: Path, tag: str) -> Path:
    """Return a sidecar path such as sample.raw.bed.gz from sample.bed.gz."""
    suffixes = path.suffixes
    if suffixes and suffixes[-1] == ".gz":
        base_suffixes = suffixes[:-1]
        gz = ".gz"
    else:
        base_suffixes = suffixes
        gz = ""
    base_name = path.name
    for suffix in suffixes:
        base_name = base_name[: -len(suffix)]
    return path.with_name(f"{base_name}{tag}{''.join(base_suffixes)}{gz}")


def _sidecar_with_suffix(path: Path, tag: str, suffix: str) -> Path:
    """Return a sidecar path with a deliberate suffix, e.g. sample.raw_nuc.bed."""
    base_name = path.name
    for existing_suffix in path.suffixes:
        base_name = base_name[: -len(existing_suffix)]
    return path.with_name(f"{base_name}{tag}{suffix}")


def _comma_ints(value: str, *, field_name: str, source: Path, line_number: int) -> list[int]:
    try:
        return [int(x) for x in value.rstrip(",").split(",") if x != ""]
    except ValueError as exc:
        raise ValueError(
            f"Invalid integer list in {field_name} at {source}:{line_number}: {value!r}"
        ) from exc


def _parse_block_bed_record(
    line: str,
    *,
    source: Path,
    line_number: int,
) -> BlockBedRecord:
    """Parse BED12/BED9+2/BED8+2 using block-list content, not only column count."""
    columns = line.rstrip("\n").split("\t")
    if len(columns) < 10:
        raise ValueError(
            f"Invalid fibertools BED row at {source}:{line_number}: "
            f"expected at least 10 columns, found {len(columns)}"
        )
    try:
        chrom_start = int(columns[1])
        chrom_end = int(columns[2])
        thick_start = int(columns[6])
        thick_end = int(columns[7])
    except ValueError as exc:
        raise ValueError(
            f"Invalid coordinate field at {source}:{line_number}: "
            "BED start/end/thickStart/thickEnd must be integers"
        ) from exc
    if chrom_start > chrom_end or thick_start > thick_end:
        raise ValueError(f"Invalid coordinate ordering at {source}:{line_number}")
    if columns[5] not in {"+", "-", "."}:
        raise ValueError(f"Invalid strand at {source}:{line_number}: {columns[5]!r}")

    if len(columns) >= 12:
        # BED12: blockCount, blockSizes, blockStarts.
        candidates: list[tuple[int | None, int, int]] = [(9, 10, 11)]
    elif len(columns) == 11:
        # BED9+2: itemRgb, blockSizes, blockStarts.
        candidates = [(None, 9, 10)]
    else:
        # BED8+2: blockSizes, blockStarts.  Do not use this as a fallback for
        # longer malformed rows because an RGB field can look like a block list.
        candidates = [(None, 8, 9)]

    errors: list[str] = []
    for count_idx, sizes_idx, starts_idx in candidates:
        try:
            sizes = _comma_ints(
                columns[sizes_idx],
                field_name="blockSizes",
                source=source,
                line_number=line_number,
            )
            starts = _comma_ints(
                columns[starts_idx],
                field_name="blockStarts",
                source=source,
                line_number=line_number,
            )
            if not sizes or not starts:
                raise ValueError("empty block list")
            if len(sizes) != len(starts):
                raise ValueError(f"{len(sizes)} sizes != {len(starts)} starts")
            if count_idx is not None and int(columns[count_idx]) != len(sizes):
                raise ValueError(
                    f"blockCount {columns[count_idx]} != {len(sizes)} parsed blocks"
                )
            span = chrom_end - chrom_start
            if any(size < 0 for size in sizes) or any(start < 0 for start in starts):
                raise ValueError("negative block size/start")
            if any(start > span for start in starts):
                raise ValueError("block start outside read span")
            return BlockBedRecord(columns, count_idx, sizes_idx, starts_idx, sizes, starts)
        except (ValueError, IndexError) as exc:
            errors.append(str(exc))

    raise ValueError(
        f"Could not identify consistent block columns at {source}:{line_number}. "
        f"Tried BED12, BED9+2, and BED8+2. Details: {'; '.join(errors)}"
    )


def normalize_fibertools_block_bed(
    input_bed: Path,
    output_bed: Path,
    *,
    expected: str = "any",
) -> dict[str, int | str]:
    """Validate fibertools BED-like output and write builder-compatible BED12.

    Args:
        input_bed: Raw ``ft extract`` output.
        output_bed: Normalized BED12 path for the builder.
        expected: ``single_base`` for 6mA-like calls, ``interval`` for MSP/nuc,
            or ``any`` when only format normalization is desired.

    Returns:
        Small QC summary useful for the extraction manifest.
    """
    if expected not in {"any", "single_base", "interval"}:
        raise ValueError("expected must be one of: any, single_base, interval")

    output_bed.parent.mkdir(parents=True, exist_ok=True)
    rows = 0
    blocks = 0
    single_base_blocks = 0
    with smart_open(input_bed) as source, smart_open(output_bed, "wt") as target:
        for line_number, line in enumerate(source, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            record = _parse_block_bed_record(
                line,
                source=input_bed,
                line_number=line_number,
            )
            inner_sizes = record.sizes[1:-1] if len(record.sizes) > 2 else record.sizes
            single_base_blocks += sum(1 for size in inner_sizes if size == 1)
            blocks += len(record.sizes)
            if expected == "single_base" and inner_sizes:
                non_single = sum(1 for size in inner_sizes if size != 1)
                if non_single / len(inner_sizes) > 0.05:
                    raise ValueError(
                        f"{input_bed}:{line_number} does not look like a single-base "
                        f"track: {non_single}/{len(inner_sizes)} internal blocks have size != 1"
                    )

            columns = record.columns
            if record.block_count_index is None:
                # Convert BED9+2 or BED8+2 to standard BED12.  BED8+2 receives a
                # placeholder itemRgb so the resulting columns are canonical.
                if record.sizes_index == 8:
                    columns = columns[:8] + ["0", str(len(record.sizes))] + columns[8:]
                else:
                    columns = columns[:9] + [str(len(record.sizes))] + columns[9:]
            target.write("\t".join(columns[:12]) + "\n")
            rows += 1

    if rows == 0:
        raise ValueError(f"No data rows found in fibertools BED output: {input_bed}")
    return {
        "input": str(input_bed),
        "output": str(output_bed),
        "rows": rows,
        "blocks": blocks,
        "single_base_blocks": single_base_blocks,
    }


def _read_all_header(path: Path) -> tuple[list[str], dict[str, int]]:
    with smart_open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
    normalized = [field[1:] if field.startswith("#") else field for field in header]
    return normalized, {field: i for i, field in enumerate(normalized)}


def _all_list(value: str) -> list[str]:
    return [x for x in value.rstrip(",").split(",") if x != ""]


def inspect_pacbio_all(
    all_tsv: Path,
    *,
    max_records: int = 1_000,
) -> dict[str, int | str | bool]:
    """Validate named ``ft extract --all`` columns and list alignment.

    This is intentionally header-driven.  If fibertools changes column order but
    keeps the field names, validation and conversion still work.  If field names
    disappear, the error tells us exactly which semantic input is missing.
    """
    header, index = _read_all_header(all_tsv)
    missing_5mc = sorted(_ALL_REQUIRED_5MC_COLUMNS - set(index))
    if missing_5mc:
        raise ValueError(
            f"{all_tsv} is missing required PacBio 5mC columns: {missing_5mc}. "
            f"Observed columns: {header}"
        )
    has_fire = _ALL_FIRE_COLUMNS <= set(index)

    records = 0
    mapped_5mc = 0
    total_5mc = 0
    mapped_m6a = 0
    total_m6a = 0
    fire_entries = 0
    positive_fire_entries = 0
    with smart_open(all_tsv) as f:
        _ = f.readline()
        for line_number, line in enumerate(f, start=2):
            if not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < len(header):
                raise ValueError(
                    f"Short row in {all_tsv}:{line_number}: "
                    f"expected {len(header)} fields, found {len(cols)}"
                )
            records += 1
            for required in ["ct", "st", "en", "fiber", "strand"]:
                if cols[index[required]] == "":
                    raise ValueError(f"Empty {required!r} field at {all_tsv}:{line_number}")
            int(cols[index["st"]])
            int(cols[index["en"]])

            ref_5mc = _all_list(cols[index["ref_5mC"]])
            qual_5mc = _all_list(cols[index["5mC_qual"]])
            if len(ref_5mc) != len(qual_5mc):
                raise ValueError(
                    f"5mC list mismatch at {all_tsv}:{line_number}: "
                    f"{len(ref_5mc)} ref_5mC != {len(qual_5mc)} 5mC_qual"
                )
            for pos, qual in zip(ref_5mc, qual_5mc, strict=True):
                pos_i = int(pos)
                qual_i = int(qual)
                if not 0 <= qual_i <= 255:
                    raise ValueError(
                        f"5mC_qual out of ML-byte range at {all_tsv}:{line_number}: {qual_i}"
                    )
                total_5mc += 1
                if pos_i >= 0:
                    mapped_5mc += 1

            if {"ref_m6a", "m6a_qual"} <= set(index):
                ref_m6a = _all_list(cols[index["ref_m6a"]])
                qual_m6a = _all_list(cols[index["m6a_qual"]])
                if len(ref_m6a) != len(qual_m6a):
                    raise ValueError(
                        f"m6A list mismatch at {all_tsv}:{line_number}: "
                        f"{len(ref_m6a)} ref_m6a != {len(qual_m6a)} m6a_qual"
                    )
                total_m6a += len(ref_m6a)
                mapped_m6a += sum(1 for pos in ref_m6a if int(pos) >= 0)

            if has_fire:
                starts = _all_list(cols[index["ref_msp_starts"]])
                lengths = _all_list(cols[index["ref_msp_lengths"]])
                fire = _all_list(cols[index["fire"]])
                if not (len(starts) == len(lengths) == len(fire)):
                    raise ValueError(
                        f"FIRE/MSP list mismatch at {all_tsv}:{line_number}: "
                        f"{len(starts)} starts, {len(lengths)} lengths, {len(fire)} fire scores"
                    )
                fire_entries += len(fire)
                positive_fire_entries += sum(1 for value in fire if float(value) > 0)

            if records >= max_records:
                break

    if records == 0:
        raise ValueError(f"No data rows found in PacBio all TSV: {all_tsv}")
    return {
        "path": str(all_tsv),
        "records_checked": records,
        "has_fire_columns": has_fire,
        "total_5mC_calls_checked": total_5mc,
        "mapped_5mC_calls_checked": mapped_5mc,
        "total_m6A_calls_checked": total_m6a,
        "mapped_m6A_calls_checked": mapped_m6a,
        "fire_entries_checked": fire_entries,
        "positive_fire_entries_checked": positive_fire_entries,
    }


def convert_pacbio_all_5mc_to_modkit(
    all_tsv: Path,
    output_tsv: Path,
    *,
    max_records_for_validation: int = 1_000,
) -> dict[str, int | str | bool]:
    """Convert PacBio ``ft extract --all`` 5mC calls to builder-compatible TSV.

    The builder already knows how to parse modkit-like 5mC rows, so this writes a
    minimal 14-column table with the same semantic fields used by the builder:
    read id, reference position, chromosome, normalized probability, and mod code.
    PacBio ML byte qualities are converted to 0-1 probabilities by ``qual / 255``.
    """
    summary = inspect_pacbio_all(all_tsv, max_records=max_records_for_validation)
    header, index = _read_all_header(all_tsv)
    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    rows = 0
    skipped_unmapped = 0
    with smart_open(all_tsv) as source, smart_open(output_tsv, "wt") as target:
        _ = source.readline()
        writer = csv.writer(target, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "read_id",
                "forward_read_position",
                "ref_position",
                "chrom",
                "mod_strand",
                "ref_strand",
                "ref_mod_strand",
                "fw_soft_clipped_start",
                "fw_soft_clipped_end",
                "read_length",
                "base_qual",
                "ref_kmer",
                "mod_qual",
                "mod_code",
            ]
        )
        for line in source:
            if not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            chrom = cols[index["ct"]]
            fiber = cols[index["fiber"]]
            strand = cols[index["strand"]]
            read_length = cols[index["fiber_length"]] if "fiber_length" in index else ""
            ref_positions = _all_list(cols[index["ref_5mC"]])
            qualities = _all_list(cols[index["5mC_qual"]])
            if len(ref_positions) != len(qualities):
                raise ValueError(
                    f"5mC list mismatch while converting {all_tsv}: "
                    f"{len(ref_positions)} positions != {len(qualities)} qualities"
                )
            for ref_pos, qual in zip(ref_positions, qualities, strict=True):
                pos = int(ref_pos)
                if pos < 0:
                    skipped_unmapped += 1
                    continue
                probability = int(qual) / 255.0
                writer.writerow(
                    [
                        fiber,
                        "",
                        pos,
                        chrom,
                        strand,
                        strand,
                        strand,
                        "",
                        "",
                        read_length,
                        "",
                        "",
                        f"{probability:.6f}",
                        "m",
                    ]
                )
                rows += 1

    return {
        **summary,
        "output": str(output_tsv),
        "converted_5mC_rows": rows,
        "skipped_unmapped_5mC_rows": skipped_unmapped,
    }


def _validate_outputs(sample: SampleConfig) -> dict[str, dict[str, str | int]]:
    outputs: dict[str, dict[str, str | int]] = {}
    for layer, path in sample.layers.items():
        if not path.exists():
            raise FileNotFoundError(f"Configured {layer} output was not created: {path}")
        size = path.stat().st_size
        if size == 0:
            raise ValueError(f"Configured {layer} output is empty: {path}")
        outputs[layer] = {"path": str(path), "size_bytes": size}
    return outputs


def _validate_sample_inputs(sample: SampleConfig) -> None:
    if sample.bam is None:
        raise ValueError(f"Sample {sample.name!r} has no BAM path in config")
    if not sample.bam.exists():
        raise FileNotFoundError(f"BAM not found for sample {sample.name!r}: {sample.bam}")
    if "5hmC" in sample.layers:
        raise ValueError(
            "PacBio extraction does not currently produce a 5hmC layer. "
            "Remove '5hmC' from this sample's layers for --platform pacbio."
        )


def extract_sample(config: Config, sample: SampleConfig) -> None:
    """Extract and normalize all configured PacBio layers for one sample."""
    _validate_sample_inputs(sample)
    settings = config.extraction
    ft = _resolve_executable(settings.ft_executable)
    assert sample.bam is not None
    command_log: list[list[str]] = []
    conversion_qc: dict[str, object] = {}

    bam_qc = {}
    samtools: str | None = None
    if settings.validate_bam:
        samtools = _resolve_executable(settings.samtools_executable)
        log.info(f"Validating BAM: {sample.bam}")
        bam_qc = validate_bam(sample.bam, "pacbio", samtools)
        log.info(
            "BAM validation passed: "
            f"{bam_qc['records_checked_for_tags']} sampled records readable"
        )

    for path in sample.layers.values():
        path.parent.mkdir(parents=True, exist_ok=True)

    if "5mC" in sample.layers:
        output = sample.layers["5mC"]
        raw_all = _insert_before_suffix(output, ".fibertools_all")
        if settings.overwrite or not output.exists():
            if settings.overwrite or not raw_all.exists():
                _run(
                    [
                        ft,
                        "extract",
                        "--all",
                        str(raw_all),
                        "--simplify",
                        str(sample.bam),
                    ],
                    command_log,
                )
            conversion_qc["5mC"] = convert_pacbio_all_5mc_to_modkit(raw_all, output)
            log.info(f"Wrote PacBio 5mC builder TSV: {output}")
        else:
            log.info(f"Skipping existing PacBio 5mC TSV: {output}")

    if "6mA" in sample.layers:
        output = sample.layers["6mA"]
        raw_bed = _insert_before_suffix(output, ".raw")
        if settings.overwrite or not output.exists():
            if settings.overwrite or not raw_bed.exists():
                _run([ft, "extract", "--m6a", str(raw_bed), str(sample.bam)], command_log)
            conversion_qc["6mA"] = normalize_fibertools_block_bed(
                raw_bed,
                output,
                expected="single_base",
            )
            log.info(f"Wrote normalized PacBio 6mA BED12: {output}")
        else:
            log.info(f"Skipping existing PacBio 6mA BED12: {output}")

    if "msp" in sample.layers:
        output = sample.layers["msp"]
        raw_bed = _insert_before_suffix(output, ".raw")
        if settings.overwrite or not output.exists():
            if settings.overwrite or not raw_bed.exists():
                _run([ft, "extract", "--msp", str(raw_bed), str(sample.bam)], command_log)
            conversion_qc["msp"] = normalize_fibertools_block_bed(
                raw_bed,
                output,
                expected="interval",
            )
            log.info(f"Wrote normalized PacBio MSP BED12: {output}")
        else:
            log.info(f"Skipping existing PacBio MSP BED12: {output}")

    nuc_csv = sample.layers["nucleosomes"]
    nuc_bed = _sidecar_with_suffix(nuc_csv, ".raw_nuc", ".bed")
    if settings.overwrite or not nuc_csv.exists():
        if settings.overwrite or not nuc_bed.exists():
            _run([ft, "extract", "--nuc", str(nuc_bed), str(sample.bam)], command_log)
        normalized_nuc_bed = _sidecar_with_suffix(nuc_csv, ".nuc_bed12", ".bed")
        conversion_qc["nucleosomes_bed"] = normalize_fibertools_block_bed(
            nuc_bed,
            normalized_nuc_bed,
            expected="interval",
        )
        rows = flatten_nucleosome_bed12(normalized_nuc_bed, nuc_csv)
        conversion_qc["nucleosomes_csv"] = {"rows": rows, "output": str(nuc_csv)}
        log.info(f"Wrote {rows:,} nucleosome rows: {nuc_csv}")
        if not settings.keep_nucleosome_bed:
            nuc_bed.unlink(missing_ok=True)
            normalized_nuc_bed.unlink(missing_ok=True)
    else:
        log.info(f"Skipping existing PacBio nucleosome CSV: {nuc_csv}")

    outputs = _validate_outputs(sample)
    if settings.write_manifest:
        executables = {"fibertools-rs": ft}
        if samtools is not None:
            executables["samtools"] = samtools
        output_dir = nuc_csv.parent
        manifest_path = write_manifest(
            output_dir,
            {
                "sample": sample.name,
                "platform": "pacbio",
                "input": {
                    "bam": str(sample.bam),
                    "bam_size_bytes": sample.bam.stat().st_size,
                },
                "bam_qc": bam_qc,
                "tools": capture_tool_versions(executables),
                "parameters": {
                    "threads": settings.threads,
                    "overwrite": settings.overwrite,
                    "methylation_threshold": config.parameters.methylation_threshold,
                },
                "commands_run": command_log,
                "conversion_qc": conversion_qc,
                "outputs": outputs,
                "output_inventory": compute_qc_metrics(output_dir),
            },
            filename=f"PACKAGE_manifest_{sample.name}.json",
        )
        log.info(f"Wrote extraction manifest: {manifest_path}")


def extract_samples(config: Config, samples: list[str] | None = None) -> None:
    """Run PacBio extraction pipeline on one or more configured samples."""
    sample_names = samples or config.sample_names
    log.info(f"PacBio extraction: samples={sample_names}")
    for name in sample_names:
        extract_sample(config, config.get_sample(name))
