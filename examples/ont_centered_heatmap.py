#!/usr/bin/env python3
"""Export center-based matrices and save paired heatmap/metaplot figures."""

from __future__ import annotations

import argparse
from pathlib import Path

from PACKAGE.analysis import export_centered_annotation_matrices
from PACKAGE.analysis.heatmap import EXTENSION_CONFIG, PROMOTER_TYPES
from PACKAGE.viz import plot_centered_annotation_heatmap, plot_centered_metaplot


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--annotation", required=True)
    parser.add_argument("--samples", nargs="+", required=True)
    parser.add_argument("--bed", type=Path, help="optional 9-column annotation BED for strand")
    parser.add_argument("--max-regions", type=int)
    parser.add_argument("--min-fibers", type=int, default=10)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    export_centered_annotation_matrices(
        args.db,
        args.outdir,
        annotations=[args.annotation],
        samples=args.samples,
        bed_path=args.bed,
        max_regions=args.max_regions,
        min_fibers=args.min_fibers,
    )

    extension_bp = EXTENSION_CONFIG[args.annotation]
    center_label = "TSS" if args.annotation in PROMOTER_TYPES else "Center"
    metaplot_inputs: dict[str, Path] = {}
    for sample in args.samples:
        prefix = f"{args.annotation}_{sample}"
        heatmap = plot_centered_annotation_heatmap(
            args.outdir / f"{prefix}_nuc.npz",
            args.outdir / f"{prefix}_met.npz",
            title=f"{args.annotation} ({sample}) centered at {center_label}",
        )
        heatmap.savefig(args.outdir / f"{prefix}_center_heatmap.pdf")
        heatmap.savefig(args.outdir / f"{prefix}_center_heatmap.png", dpi=200)
        metaplot_inputs[sample] = args.outdir / f"{prefix}_metaplot.csv"

    metaplot = plot_centered_metaplot(
        metaplot_inputs,
        annotation=args.annotation,
        extension_bp=extension_bp,
        center_label=center_label,
    )
    metaplot.savefig(args.outdir / f"{args.annotation}_center_metaplot.pdf")
    metaplot.savefig(args.outdir / f"{args.annotation}_center_metaplot.png", dpi=200)
    print(f"Wrote centered heatmap outputs to {args.outdir}")


if __name__ == "__main__":
    main()
