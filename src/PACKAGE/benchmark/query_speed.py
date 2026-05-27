"""Query speed benchmark.

Measures wall-clock time for region queries as a function of region size, across a
range of sizes (e.g. 1 kb, 10 kb, 100 kb, 1 Mb). Compares against a naive baseline
(loading the entire BED into pandas + intersection) when available.

STATUS: Stub.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any


def run_query_speed_benchmark(
    db: Any,
    output_dir: Path,
    sample: str | None = None,
    region_sizes_kb: list[int] | None = None,
    n_replicates: int = 10,
) -> None:
    """Run the query-speed benchmark and write outputs.

    Outputs:
        - ``query_speed_results.csv``: per-region-size timings
        - ``query_speed.pdf``: figure (paper Fig 2a or extended data)
    """
    raise NotImplementedError
