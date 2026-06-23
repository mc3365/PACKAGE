#!/usr/bin/env python3
"""Export global per-fiber feature fractions and plot ECDF panels."""

from __future__ import annotations

import argparse
from pathlib import Path

from PACKAGE.analysis import export_global_feature_fractions
from PACKAGE.viz import plot_feature_ecdfs


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--samples", nargs="+")
    parser.add_argument("--max-fibers-per-chrom", type=int)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    feature_csv = args.outdir / "global_feature_fractions.csv"
    export_global_feature_fractions(
        args.db,
        feature_csv,
        samples=args.samples,
        max_fibers_per_chrom=args.max_fibers_per_chrom,
    )
    coverage_fig, base_fig = plot_feature_ecdfs(feature_csv)
    coverage_fig.savefig(args.outdir / "ecdf_coverage_fraction.pdf")
    coverage_fig.savefig(args.outdir / "ecdf_coverage_fraction.png", dpi=200)
    base_fig.savefig(args.outdir / "ecdf_base_specific_fraction.pdf")
    base_fig.savefig(args.outdir / "ecdf_base_specific_fraction.png", dpi=200)
    print(f"Wrote ECDF outputs to {args.outdir}")


if __name__ == "__main__":
    main()
