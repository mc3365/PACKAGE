"""Command-line interface for PACKAGE.

Exposed after `pip install` as the `PACKAGE` command. Subcommands:

    PACKAGE extract    Extract modifications from BAM files (ONT or PacBio)
    PACKAGE build      Build an HDF5 database from extracted intermediate files
    PACKAGE query      Query an existing HDF5 database
    PACKAGE info       Print database summary

This module is intentionally thin — it parses arguments and dispatches to library
functions. All logic lives in the library modules. Scripting users should bypass the
CLI and import the library directly.
"""

from __future__ import annotations

from pathlib import Path

import click

from PACKAGE import __version__


@click.group()
@click.version_option(version=__version__, prog_name="PACKAGE")
def cli() -> None:
    """PACKAGE: single-molecule epigenomic analysis for long-read Fiber-seq data."""


# ---------------------------------------------------------------------------
# `PACKAGE extract` — extraction from BAM
# ---------------------------------------------------------------------------
@cli.command()
@click.option(
    "--config",
    "config_path",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    required=True,
    help="Path to YAML config file describing inputs and outputs.",
)
@click.option(
    "--platform",
    type=click.Choice(["ont", "pacbio"], case_sensitive=False),
    required=True,
    help="Sequencing platform of the input BAM files.",
)
@click.option(
    "--samples",
    multiple=True,
    help="Sample names to extract (default: all in config). Repeat flag for multiple.",
)
def extract(config_path: Path, platform: str, samples: tuple[str, ...]) -> None:
    """Extract modification calls and accessibility features from BAM files.

    This wraps modkit (5mC/5hmC) + fibertools-rs ft extract (m6A, MSP, nucleosomes) for
    ONT data, or pb-CpG-tools + fibertools-rs for PacBio data. Outputs intermediate
    files that can then be passed to `PACKAGE build`.
    """
    from PACKAGE.config import load_config

    cfg = load_config(config_path)
    sample_list = list(samples) if samples else cfg.sample_names
    click.echo(f"Extracting for platform={platform}, samples={sample_list}")

    if platform.lower() == "ont":
        from PACKAGE.extract.ont import extract_samples

        extract_samples(cfg, samples=sample_list)
    else:
        from PACKAGE.extract.pacbio import extract_samples

        extract_samples(cfg, samples=sample_list)
    click.secho("✓ Extraction complete.", fg="green")


# ---------------------------------------------------------------------------
# `PACKAGE build` — HDF5 database build
# ---------------------------------------------------------------------------
@cli.command()
@click.option(
    "--config",
    "config_path",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    required=True,
    help="Path to YAML config file describing inputs and outputs.",
)
@click.option(
    "--samples",
    multiple=True,
    help="Sample names to build (default: all in config). Repeat flag for multiple.",
)
def build(config_path: Path, samples: tuple[str, ...]) -> None:
    """Build an HDF5 database from extracted intermediate files.

    Expects that `PACKAGE extract` (or an equivalent manual extraction) has produced
    the intermediate files described in the config.
    """
    from PACKAGE.config import load_config
    from PACKAGE.db.builder import build_database

    cfg = load_config(config_path)
    sample_list = list(samples) if samples else cfg.sample_names
    click.echo(f"Building database {cfg.output_path} for samples={sample_list}")
    build_database(cfg, samples=sample_list)
    click.secho("✓ Build complete.", fg="green")


# ---------------------------------------------------------------------------
# `PACKAGE query` — query database
# ---------------------------------------------------------------------------
@cli.command()
@click.option(
    "--db",
    "db_path",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    required=True,
    help="Path to the HDF5 database.",
)
@click.option("--region", help="Region as chrom:start-end (e.g. chr1:1000-2000).")
@click.option("--sample", default=None, help="Restrict to one sample.")
@click.option(
    "--out",
    "out_path",
    type=click.Path(path_type=Path),
    default=None,
    help="Optional output CSV path. If omitted, prints to stdout.",
)
def query(
    db_path: Path,
    region: str | None,
    sample: str | None,
    out_path: Path | None,
) -> None:
    """Query a PACKAGE database for fibers and their modifications in a region."""
    from PACKAGE.db.database import FiberDatabase

    if not region:
        click.secho("Error: --region is required (e.g. chr1:1000-2000).", fg="red")
        raise SystemExit(1)
    chrom, span = region.split(":")
    start_str, end_str = span.split("-")
    start, end = int(start_str), int(end_str)

    with FiberDatabase(db_path) as db:
        fibers = db.get_fibers_at(chrom, start, end, sample=sample)
    click.echo(f"Found {len(fibers)} fibers overlapping {region}")
    if out_path:
        out_path.write_text("\n".join(fibers))
        click.echo(f"Wrote to {out_path}")


# ---------------------------------------------------------------------------
# `PACKAGE info` — quick database summary
# ---------------------------------------------------------------------------
@cli.command()
@click.argument(
    "db_path",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
def info(db_path: Path) -> None:
    """Print a summary of a PACKAGE HDF5 database."""
    from PACKAGE.db.database import FiberDatabase

    with FiberDatabase(db_path) as db:
        summary = db.get_summary()
    for key, value in summary.items():
        click.echo(f"  {key}: {value}")


if __name__ == "__main__":
    cli()
