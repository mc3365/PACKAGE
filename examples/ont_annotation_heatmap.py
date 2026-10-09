#!/usr/bin/env python3
"""Export annotation-centered matrices and save heatmap/metaplot figures."""

from __future__ import annotations

import argparse
from pathlib import Path

from mei_fiber.analysis import export_annotation_matrices
from mei_fiber.viz import plot_annotation_heatmap, plot_annotation_metaplot


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--annotation", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--max-regions", type=int)
    parser.add_argument("--min-fibers", type=int, default=10)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    export_annotation_matrices(
        args.db,
        args.outdir,
        annotations=[args.annotation],
        samples=[args.sample],
        max_regions=args.max_regions,
        min_fibers=args.min_fibers,
    )

    prefix = f"{args.annotation}_{args.sample}"
    heatmap = plot_annotation_heatmap(
        args.outdir / f"{prefix}_nuc_body.npz",
        args.outdir / f"{prefix}_met_body.npz",
        title=f"{args.sample} {args.annotation}",
    )
    metaplot = plot_annotation_metaplot(
        args.outdir / f"{prefix}_nuc_metaplot.npz",
        args.outdir / f"{prefix}_met_metaplot.npz",
        title=f"{args.sample} {args.annotation}",
    )
    heatmap.savefig(args.outdir / f"{prefix}_heatmap.pdf")
    heatmap.savefig(args.outdir / f"{prefix}_heatmap.png", dpi=200)
    metaplot.savefig(args.outdir / f"{prefix}_metaplot.pdf")
    metaplot.savefig(args.outdir / f"{prefix}_metaplot.png", dpi=200)
    print(f"Wrote annotation matrix outputs to {args.outdir}")


if __name__ == "__main__":
    main()
