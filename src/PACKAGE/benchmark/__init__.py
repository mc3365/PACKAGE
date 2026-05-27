"""Benchmark suite for the methods paper.

The benchmarks here reproduce the paper's performance and validation figures. Each
submodule writes its outputs (numeric tables + matplotlib figures) to a configurable
output directory so the paper figures can be regenerated end-to-end from a built
database.

Submodules:
    query_speed     Query time vs region size (paper Fig 2a)
    memory_profile  Memory footprint vs database size
    build_time      Database build time profiling
    cross_platform  ONT vs PacBio concordance metrics (paper Fig 2b)

All benchmarks take an open ``FiberDatabase`` and write outputs to ``--outdir``.
"""

__all__ = []  # nothing re-exported; users import submodules directly
