"""Spatial index for fast region-based fiber queries.

Builds an ``IntervalTree`` per (sample, chromosome) keyed on fiber spans. The tree
stores uint32 fiber IDs as payloads; the caller resolves them to UUIDs via the
per-chromosome ``fiber_id_table``.

Monster fibers (span > ``FiberDatabase.MAX_FIBER_SPAN``) are excluded at index-build
time. This means the index always returns "clean" overlap sets without needing
post-filtering. The threshold passed into ``build_for_sample`` should match the
``FiberDatabase`` class constant so the spatial index and the array-scan fallback
agree.

The index is persisted as a sidecar pickle (``<db>.index.pkl``) alongside the HDF5
file. Loading takes O(seconds); building takes O(minutes) on a full mouse-genome
dataset, which is why it's cached.

This is a refactor of the spatial-index code from the V8 ``FiberDatabase`` class,
moved to its own module so the database reader can stay focused on querying.
"""

from __future__ import annotations

import pickle
from collections.abc import Callable
from pathlib import Path
from typing import Any

import h5py
from intervaltree import IntervalTree

from PACKAGE.db import schema
from PACKAGE.utils import get_logger

log = get_logger(__name__)


def build_for_sample(
    db: h5py.File,
    sample: str,
    get_id_maps: Callable[[str, str], tuple[Any, Any]],
    max_fiber_span: int = 100_000,
) -> dict[str, IntervalTree]:
    """Build an IntervalTree per chromosome for one sample.

    Args:
        db: Open ``h5py.File`` handle.
        sample: Sample name (top-level HDF5 group).
        get_id_maps: Callable ``(sample, chrom) -> (str_to_int, int_to_str_table)``
            from the FiberDatabase that owns this build. Passed in (rather than
            recomputed here) so the caches are shared.
        max_fiber_span: Fibers with span > this are excluded as likely artifacts.
            Default matches ``FiberDatabase.MAX_FIBER_SPAN``. Pass a different
            value only if you've subclassed the database and changed the threshold.

    Returns:
        ``{chrom: IntervalTree}`` mapping. Each tree's intervals carry the uint32
        fiber ID as the ``data`` payload. Monster fibers are excluded.
    """
    # We only know which chromosomes exist by inspecting the sample group.
    chroms = sorted(k for k in db[sample].keys() if k.startswith("chr"))
    result: dict[str, IntervalTree] = {}
    for i, chrom in enumerate(chroms, 1):
        meta_path = schema.fiber_metadata_path(sample, chrom)
        fids = db[f"{meta_path}/fiber_int_ids"][:]
        starts = db[f"{meta_path}/starts"][:]
        ends = db[f"{meta_path}/ends"][:]
        # Filter monster fibers (artifact entries with implausibly large spans).
        spans = ends - starts
        valid = spans <= max_fiber_span
        n_excluded = int((~valid).sum())
        if n_excluded > 0:
            fids, starts, ends = fids[valid], starts[valid], ends[valid]
        # Touch the id_maps so any caching the caller wants happens up front.
        get_id_maps(sample, chrom)
        tree = IntervalTree()
        for int_id, s, e in zip(fids, starts, ends, strict=False):
            if int(e) > int(s):  # IntervalTree requires non-empty intervals
                tree.addi(int(s), int(e), int(int_id))
        result[chrom] = tree
        suffix = f" (excluded {n_excluded:,} monster fibers)" if n_excluded else ""
        log.info(f"  [{sample}] {chrom} ({i}/{len(chroms)}): {len(tree)} fibers{suffix}")
    return result


def save_index(index: dict[str, dict[str, IntervalTree]], path: Path) -> None:
    """Persist a spatial index to disk as a pickle sidecar."""
    with path.open("wb") as f:
        pickle.dump(index, f)


def load_index(path: Path) -> dict[str, dict[str, IntervalTree]]:
    """Load a previously-saved spatial index from disk."""
    with path.open("rb") as f:
        return pickle.load(f)
