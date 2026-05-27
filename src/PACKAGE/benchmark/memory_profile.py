"""Memory footprint benchmark.

Profiles peak RSS during representative operations (open database, build spatial
index, run a large region query) as a function of database size.

STATUS: Stub.
"""

from __future__ import annotations

from pathlib import Path


def run_memory_benchmark(db_paths: list[Path], output_dir: Path) -> None:
    """Profile memory across multiple database sizes."""
    raise NotImplementedError
