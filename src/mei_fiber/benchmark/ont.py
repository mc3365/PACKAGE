"""Command-line entry point for the ONT storage and query benchmarks."""

from __future__ import annotations

import argparse
from pathlib import Path

from mei_fiber.db import FiberDatabase

from .query_speed import run_query_speed_benchmark
from .storage import run_storage_benchmark


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db", type=Path, required=True, help="MEI-Fiber HDF5 database")
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--sample", help="sample to benchmark; defaults to the first")
    parser.add_argument("--config", type=Path, help="optional YAML for intermediate sizes")
    parser.add_argument("--region-sizes-kb", type=int, nargs="+", default=[1, 10, 100, 1000])
    parser.add_argument("--n-regions", type=int, default=25)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--seed", type=int, default=3365)
    parser.add_argument(
        "--query-modes",
        nargs="+",
        choices=["array_scan", "spatial_index"],
        default=["array_scan", "spatial_index"],
        help=(
            "query implementations to benchmark; use 'spatial_index' alone for "
            "large production databases"
        ),
    )
    parser.add_argument(
        "--build-index",
        action="store_true",
        help="build and save the spatial index before benchmarking",
    )
    args = parser.parse_args()

    run_storage_benchmark(args.db, args.outdir, config_path=args.config)
    with FiberDatabase(args.db, build_index=args.build_index) as db:
        run_query_speed_benchmark(
            db,
            args.outdir,
            sample=args.sample,
            region_sizes_kb=args.region_sizes_kb,
            n_regions=args.n_regions,
            repeats=args.repeats,
            seed=args.seed,
            query_modes=args.query_modes,
        )
    print(f"Wrote ONT benchmark outputs to {args.outdir}")


if __name__ == "__main__":
    main()
