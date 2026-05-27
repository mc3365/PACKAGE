"""Spatial index for fast region-based fiber queries.

Builds an ``IntervalTree`` of fiber spans per (sample, chromosome). The tree stores
integer fiber IDs as payloads; the caller resolves them to UUIDs via the per-chromosome
``fiber_id_table``.

The index is persisted as a sidecar pickle (``<db>.index.pkl``) alongside the HDF5 file.
Loading takes ~1-2s; building takes ~minutes on a full mouse genome dataset.

STATUS: Stub. TODO(refactor): port from fiber_database_v8.FiberDatabase._build_spatial_index.
"""

from __future__ import annotations

import pickle
from pathlib import Path
from typing import Any


def build_spatial_index(db: Any, sample: str) -> dict[str, Any]:
    """Build an IntervalTree per chromosome for one sample.

    Args:
        db: Open ``h5py.File`` handle.
        sample: Sample name (top-level HDF5 group).

    Returns:
        ``{chrom: IntervalTree}`` mapping.
    """
    raise NotImplementedError(
        "Port from fiber_database_v8.FiberDatabase._build_spatial_index"
    )


def save_spatial_index(index: dict, path: Path) -> None:
    """Persist a spatial index to disk as pickle."""
    with path.open("wb") as f:
        pickle.dump(index, f)


def load_spatial_index(path: Path) -> dict:
    """Load a previously saved spatial index from disk."""
    with path.open("rb") as f:
        return pickle.load(f)
