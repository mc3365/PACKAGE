#!/usr/bin/env python3
"""Create a single-molecule regional plot from a MEI-Fiber ONT database."""

from __future__ import annotations

import argparse
from pathlib import Path

from mei_fiber.db import FiberDatabase
from mei_fiber.viz import single_molecule_heatmap


def parse_region(value: str) -> tuple[str, int, int]:
    """Parse ``chrom:start-end`` into its components."""
    try:
        chrom, span = value.split(":", maxsplit=1)
        start, end = (int(part.replace(",", "")) for part in span.split("-", maxsplit=1))
    except ValueError as exc:
        raise argparse.ArgumentTypeError("region must be chrom:start-end") from exc
    if start < 0 or end <= start:
        raise argparse.ArgumentTypeError("region must satisfy 0 <= start < end")
    return chrom, start, end


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db", required=True, type=Path, help="MEI-Fiber HDF5 database")
    parser.add_argument("--region", required=True, type=parse_region, help="chrom:start-end")
    parser.add_argument("--sample", help="sample name; defaults to the first sample")
    parser.add_argument("--layer", choices=("5mC", "5hmC", "6mA"), default="5mC")
    parser.add_argument(
        "--tracks",
        choices=("full", "modification", "chromatin"),
        default="full",
        help=(
            "track set to draw: full shows chromatin plus the selected modification; "
            "modification shows only the selected modification; chromatin shows "
            "nucleosomes and MSPs only"
        ),
    )
    parser.add_argument("--max-fibers", type=int, default=100)
    parser.add_argument(
        "--hide-unmodified",
        action="store_true",
        help="hide unmodified 5mC/5hmC calls for a cleaner modification-focused view",
    )
    parser.add_argument("--out", type=Path, default=Path("ont_region.png"))
    args = parser.parse_args()

    chrom, start, end = args.region
    with FiberDatabase(args.db) as db:
        fig = single_molecule_heatmap(
            db,
            chrom,
            start,
            end,
            sample=args.sample,
            layer=args.layer,
            max_fibers=args.max_fibers,
            show_unmodified=not args.hide_unmodified,
            tracks=args.tracks,
        )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=200, bbox_inches="tight")
    print(f"Wrote {args.out}")


if __name__ == "__main__":
    main()
