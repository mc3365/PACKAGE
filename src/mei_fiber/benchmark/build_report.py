"""Create a machine-readable report for an instrumented full database build."""

from __future__ import annotations

import argparse
import json
import os
import platform
import socket
from pathlib import Path
from typing import Any

import h5py

from mei_fiber import __version__
from mei_fiber.config import load_config
from mei_fiber.db import schema


def parse_gnu_time(path: Path) -> dict[str, str | int | float]:
    """Parse stable fields from GNU ``/usr/bin/time -v`` output."""
    fields: dict[str, str | int | float] = {}
    for line in Path(path).read_text().splitlines():
        if ": " not in line:
            continue
        key, value = (part.strip() for part in line.split(": ", maxsplit=1))
        normalized = key.lower().replace(" ", "_").replace("%", "percent")
        try:
            fields[normalized] = int(value)
        except ValueError:
            try:
                fields[normalized] = float(value)
            except ValueError:
                fields[normalized] = value
    return fields


def write_build_report(
    config_path: Path,
    resource_log: Path,
    output_path: Path,
    started_at: str,
    finished_at: str,
    elapsed_seconds: float,
    git_commit: str,
) -> dict[str, Any]:
    """Write build metadata and resource measurements to JSON."""
    config = load_config(config_path)
    db_path = config.output_path
    if not db_path.exists():
        raise FileNotFoundError(f"Completed database not found: {db_path}")

    with h5py.File(db_path, "r") as db:
        samples = [name for name in db if name not in schema.RESERVED_TOP_LEVEL]
        sample_fibers: dict[str, int] = {}
        for sample in samples:
            if "n_fibers" in db[sample].attrs:
                sample_fibers[sample] = int(db[sample].attrs["n_fibers"])
                continue
            chroms = [name for name in db[sample] if name.startswith("chr")]
            sample_fibers[sample] = sum(
                len(db[f"{schema.fiber_metadata_path(sample, chrom)}/fiber_int_ids"])
                for chrom in chroms
            )

    slurm_keys = [
        "SLURM_JOB_ID",
        "SLURM_JOB_NAME",
        "SLURM_JOB_PARTITION",
        "SLURM_CPUS_PER_TASK",
        "SLURM_MEM_PER_NODE",
        "SLURM_JOB_NODELIST",
    ]
    report = {
        "mei_fiber_version": __version__,
        "git_commit": git_commit,
        "config": str(Path(config_path).resolve()),
        "database": str(db_path.resolve()),
        "database_size_bytes": db_path.stat().st_size,
        "samples": sample_fibers,
        "started_at": started_at,
        "finished_at": finished_at,
        "elapsed_seconds": elapsed_seconds,
        "host": socket.gethostname(),
        "platform": platform.platform(),
        "python": platform.python_version(),
        "slurm": {key: os.environ[key] for key in slurm_keys if key in os.environ},
        "gnu_time": parse_gnu_time(resource_log),
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--resource-log", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--started-at", required=True)
    parser.add_argument("--finished-at", required=True)
    parser.add_argument("--elapsed-seconds", type=float, required=True)
    parser.add_argument("--git-commit", required=True)
    args = parser.parse_args()
    write_build_report(
        args.config,
        args.resource_log,
        args.out,
        args.started_at,
        args.finished_at,
        args.elapsed_seconds,
        args.git_commit,
    )
    print(f"Wrote build report to {args.out}")


if __name__ == "__main__":
    main()
