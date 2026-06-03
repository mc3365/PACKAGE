"""Read-only query interface for PACKAGE HDF5 databases.

This is the **library face** of the storage layer. It is a refactor of the existing
``fiber_database_v8.py`` with the following changes from V8:

  - All HDF5 paths come from :mod:`PACKAGE.db.schema` instead of inline string literals
  - The spatial-index code is split into :mod:`PACKAGE.db.spatial_index`
  - FIRE-specific methods removed: FIRE peaks are stored as regular annotations in the
    master annotation file (consistent with the package's design that no annotation
    category gets special schema treatment)
  - Backward-compatible reading of databases built with the V8 ``master_v3/`` path
  - One small correctness fix: the array-scan fallback in :meth:`get_fibers_at` now
    uses strict half-open overlap semantics, matching the IntervalTree code path

The user-facing API is otherwise identical to V8: string fiber IDs in, string fiber
IDs out. Integer ID mapping is internal.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import h5py
import numpy as np
import pandas as pd

from PACKAGE.db import schema
from PACKAGE.utils import get_logger

log = get_logger(__name__)


# Path used by V8 builders for the master annotation set. We try the new path first,
# then fall back to this for backward compatibility with V8 databases.
_V8_ANNOTATION_FEATURES_PATH = "annotations/master_v3/features"


class FiberDatabase:
    """Read-only handle to a PACKAGE HDF5 database.

    Use as a context manager to ensure the underlying HDF5 file is closed:

        with FiberDatabase("fiber_database.h5") as db:
            df = db.query_annotation_fast("CGI", sample="d0")

    Parameters
    ----------
    db_path
        Path to the HDF5 database file.
    build_index
        If True, build the spatial index on open and persist it to a sidecar
        ``.index.pkl`` file. If False, attempt to load an existing sidecar.
    """

    def __init__(self, db_path: str | Path, build_index: bool = False) -> None:
        self.db_path = Path(db_path)
        if not self.db_path.exists():
            raise FileNotFoundError(f"Database not found: {self.db_path}")

        self.db = h5py.File(self.db_path, "r")
        self.samples: list[str] = [
            k for k in self.db.keys() if k not in schema.RESERVED_TOP_LEVEL
        ]

        # Read metadata. Use both new (schema.META_*) and legacy V8 attribute names so
        # this reader works against V8 databases without modification.
        if schema.metadata_path() in self.db:
            meta = self.db[schema.metadata_path()]
            self.version: str = meta.attrs.get(
                schema.META_SCHEMA_VERSION,
                meta.attrs.get("version", "unknown"),
            )
            self.genome_version: str = meta.attrs.get(
                schema.META_GENOME_VERSION, "unknown"
            )
            self.methylation_threshold: float = float(
                meta.attrs.get(schema.META_METHYLATION_THRESHOLD, 0.5)
            )
        else:
            self.version = "unknown"
            self.genome_version = "unknown"
            self.methylation_threshold = 0.5

        # Caches populated lazily as queries arrive.
        self._chromosomes: dict[str, list[str]] = {}
        self._spatial_index: dict[str, dict[str, Any]] = {}
        self._fiber_lookup_cache: dict[str, dict[bytes, bytes]] = {}
        # (sample, chrom) -> {fiber_id_bytes: uint32}
        self._int_map_cache: dict[tuple[str, str], dict[bytes, np.uint32]] = {}
        # (sample, chrom) -> np.array of UUID bytestrings (parallel to int IDs)
        self._str_map_cache: dict[tuple[str, str], np.ndarray] = {}

        # Spatial index lifecycle: either build, load, or warn.
        index_path = self.db_path.with_suffix(".index.pkl")
        if build_index:
            log.info("Building spatial index...")
            from PACKAGE.db.spatial_index import build_for_sample, save_index

            for sample in self.samples:
                self._spatial_index[sample] = build_for_sample(
                    self.db, sample, self._get_id_maps
                )
            save_index(self._spatial_index, index_path)
            log.info(f"Spatial index saved to {index_path}")
        elif index_path.exists():
            from PACKAGE.db.spatial_index import load_index

            self._spatial_index = load_index(index_path)
            log.info("Spatial index loaded")
        else:
            log.warning(
                "No spatial index found. Region queries will use array scan (slower). "
                "Run with build_index=True to build and cache the index."
            )

        # Log a one-line summary of what's in this database.
        info = ", ".join(f"{s}({self._count_fibers(s):,})" for s in self.samples)
        log.info(f"FiberDatabase opened: {self.db_path.name}, samples={info}")

    # ------------------------------------------------------------------
    # Context manager protocol
    # ------------------------------------------------------------------
    def __enter__(self) -> FiberDatabase:
        return self

    def __exit__(self, *exc: object) -> None:
        self.close()

    def close(self) -> None:
        """Close the underlying HDF5 file."""
        if self.db is not None:
            self.db.close()

    # ------------------------------------------------------------------
    # Internal: integer <-> string fiber ID mapping
    # ------------------------------------------------------------------
    # Per-chromosome fiber IDs are stored as uint32 in data arrays. The fiber_id_table
    # array (one per chromosome) translates uint32 -> UUID bytestring. We cache both
    # directions on first access.

    def _get_id_maps(
        self, sample: str, chrom: str
    ) -> tuple[dict[bytes, np.uint32], np.ndarray]:
        """Return (str_to_int, int_to_str_table) for one (sample, chromosome).

        The str_to_int dict maps UUID bytestrings to uint32. The int_to_str_table is
        the raw fiber_id_table array, indexed by uint32 fiber ID.
        """
        key = (sample, chrom)
        if key not in self._int_map_cache:
            table = self.db[schema.fiber_id_table_path(sample, chrom)][:]
            str_to_int = {fid: np.uint32(i) for i, fid in enumerate(table)}
            self._int_map_cache[key] = str_to_int
            self._str_map_cache[key] = table
        return self._int_map_cache[key], self._str_map_cache[key]

    def _str_to_int(
        self, fiber_id: str | bytes, sample: str, chrom: str
    ) -> np.uint32 | None:
        """Convert a UUID string (or bytes) to its per-chromosome uint32 ID.

        Returns None if the fiber is not present on this chromosome.
        """
        s2i, _ = self._get_id_maps(sample, chrom)
        encoded = fiber_id.encode() if isinstance(fiber_id, str) else fiber_id
        return s2i.get(encoded)

    def _int_to_str(self, int_id: np.uint32 | int, sample: str, chrom: str) -> str:
        """Convert a uint32 fiber ID back to its UUID string."""
        _, i2s = self._get_id_maps(sample, chrom)
        return i2s[int_id].decode()

    # ------------------------------------------------------------------
    # Sample / chromosome / annotation enumeration
    # ------------------------------------------------------------------
    def get_chromosomes(self, sample: str | None = None) -> list[str]:
        """List chromosome names for a sample (only those starting with 'chr')."""
        sample = sample or self.samples[0]
        if sample not in self._chromosomes:
            self._chromosomes[sample] = sorted(
                k for k in self.db[sample].keys() if k.startswith("chr")
            )
        return self._chromosomes[sample]

    def _annotation_features_root(self) -> str | None:
        """Return the HDF5 path to the annotation features group, or None if missing.

        Tries the new schema path first, then the V8 legacy path. This is what makes
        the reader compatible with databases built before the schema rename.
        """
        new_path = schema.annotation_features_path()
        if new_path in self.db:
            return new_path
        if _V8_ANNOTATION_FEATURES_PATH in self.db:
            return _V8_ANNOTATION_FEATURES_PATH
        return None

    def list_annotations(self) -> list[str]:
        """List available annotation type names (the prefixes the user grouped by)."""
        root = self._annotation_features_root()
        if root is None:
            return []
        return sorted(self.db[root].keys())

    def _count_fibers(self, sample: str) -> int:
        """Total number of fibers in one sample (across all chromosomes)."""
        path = schema.fiber_lookup_path(sample) + "/fiber_ids"
        return len(self.db[path]) if path in self.db else 0

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------
    def get_summary(self) -> dict[str, Any]:
        """Return a dict describing the database contents."""
        summary: dict[str, Any] = {
            "version": self.version,
            "genome": self.genome_version,
            "samples": self.samples,
            "n_annotations": len(self.list_annotations()),
            "has_spatial_index": bool(self._spatial_index),
        }
        for sample in self.samples:
            summary[f"n_fibers_{sample}"] = self._count_fibers(sample)
        return summary

    # ------------------------------------------------------------------
    # Query methods (Batch 2 / Batch 3 — to be added in subsequent commits)
    # ------------------------------------------------------------------
    def get_fibers_at(
        self, chrom: str, start: int, end: int, sample: str | None = None
    ) -> list[str]:
        """Return fiber IDs (as strings) overlapping ``[start, end)`` on ``chrom``.

        TODO: implement in Batch 2.
        """
        raise NotImplementedError("Will be added in Batch 2")

    def get_nucleosomes(
        self, fiber_id: str, chrom: str, sample: str | None = None
    ) -> dict[str, Any]:
        """TODO: implement in Batch 2."""
        raise NotImplementedError("Will be added in Batch 2")

    def get_methylation(
        self,
        fiber_id: str,
        chrom: str,
        mod_type: str = "5mC",
        sample: str | None = None,
    ) -> dict[str, Any]:
        """TODO: implement in Batch 2."""
        raise NotImplementedError("Will be added in Batch 2")

    def get_msp(
        self, fiber_id: str, chrom: str, sample: str | None = None
    ) -> dict[str, Any]:
        """TODO: implement in Batch 2."""
        raise NotImplementedError("Will be added in Batch 2")

    def query_annotation_fast(
        self,
        annotation_name: str,
        sample: str | None = None,
        feature_types: list[str] | None = None,
        max_regions: int | None = None,
    ) -> pd.DataFrame:
        """TODO: implement in Batch 3."""
        raise NotImplementedError("Will be added in Batch 3")
