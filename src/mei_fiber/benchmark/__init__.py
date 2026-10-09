"""Benchmark suite for the methods paper.

The benchmarks here reproduce the paper's performance and validation figures. Each
submodule writes its outputs (numeric tables + matplotlib figures) to a configurable
output directory so the paper figures can be regenerated end-to-end from a built
database.

Implemented submodules:
    ont          ONT storage and query benchmark command-line entry point
    query_speed  Warm-cache query time by region size and query mode
    storage      HDF5, intermediate-file, and record-count summaries
    coaccessibility  Legacy-versus-MEI-Fiber PacBio co-accessibility validation

Cross-platform concordance and full-build resource profiling remain later phases.
"""

__all__ = []  # nothing re-exported; users import submodules directly
