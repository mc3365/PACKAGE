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

    # ==================================================================
    # COORDINATE QUERIES (Batch 2)
    # ==================================================================

    def get_fibers_at(
        self,
        chrom: str,
        start: int,
        end: int,
        sample: str | None = None,
    ) -> list[str]:
        """Return fiber IDs (as strings) overlapping ``[start, end)`` on ``chrom``.

        Uses the spatial index if available, otherwise falls back to a vectorized
        array scan over the fiber_metadata starts/ends arrays. Both code paths use
        strict half-open overlap semantics: a fiber overlaps the region iff its
        interval ``[fiber_start, fiber_end)`` and the query interval ``[start, end)``
        have nonempty intersection.

        Args:
            chrom: Chromosome name (e.g. "chr1").
            start: Region start (0-based, inclusive).
            end: Region end (0-based, exclusive).
            sample: Sample to query. Defaults to the first sample in the database.

        Returns:
            List of fiber UUID strings (may be empty). Order is not guaranteed.
        """
        sample = sample or self.samples[0]

        # Fast path: use the spatial index when present.
        if sample in self._spatial_index and chrom in self._spatial_index[sample]:
            int_ids = [
                iv.data
                for iv in self._spatial_index[sample][chrom].overlap(start, end)
            ]
            _, i2s = self._get_id_maps(sample, chrom)
            return [i2s[iid].decode() for iid in int_ids]

        # Fallback: array scan. Returns [] if the chromosome's metadata is absent.
        meta_path = schema.fiber_metadata_path(sample, chrom)
        starts_path = f"{meta_path}/starts"
        if starts_path not in self.db:
            return []
        starts = self.db[starts_path][:]
        ends = self.db[f"{meta_path}/ends"][:]
        fids = self.db[f"{meta_path}/fiber_int_ids"][:]
        # Strict half-open overlap: fiber.end > region.start AND fiber.start < region.end.
        # (V8 used >= / <= which differs on boundary-touching cases; we match the
        #  IntervalTree semantics here so both code paths agree.)
        mask = (ends > start) & (starts < end)
        _, i2s = self._get_id_maps(sample, chrom)
        return [i2s[iid].decode() for iid in fids[mask]]

    # ==================================================================
    # PER-FIBER DATA ACCESS (Batch 2)
    # ==================================================================

    def _get_fiber_slice(
        self, fiber_id: str, chrom: str, sample: str, layer: str
    ) -> tuple[int, int] | None:
        """Return ``(start, end)`` row indices for one fiber in one layer's arrays.

        Uses the per-layer slice index stored at ``<sample>/<chrom>/_indices/<layer>_slices``.
        Returns None if any of:

          - the layer is not stored for this chromosome
          - the fiber is not present on this chromosome
          - the fiber has no rows in this layer (e.g. no nucleosomes called)
        """
        slice_path = schema.slice_index_path(sample, chrom, layer)
        if slice_path not in self.db:
            return None
        int_id = self._str_to_int(fiber_id, sample, chrom)
        if int_id is None:
            return None
        idx = self.db[slice_path][:]
        mask = idx["fiber_int_id"] == int_id
        if not mask.any():
            return None
        row = idx[mask][0]
        return (int(row["start"]), int(row["end"]))

    def get_nucleosomes(
        self, fiber_id: str, chrom: str, sample: str | None = None
    ) -> dict[str, np.ndarray]:
        """Return one fiber's nucleosomes on one chromosome.

        Returns:
            Dict with keys 'starts', 'ends', 'widths', 'linkers' — each a numpy
            array of the same length. Empty dict if the fiber is not present or
            has no nucleosomes called.
        """
        sample = sample or self.samples[0]
        sl = self._get_fiber_slice(fiber_id, chrom, sample, "nucleosomes")
        if sl is None:
            return {}
        s, e = sl
        base = schema.layer_path(sample, chrom, "nucleosomes")
        return {
            "starts":  self.db[f"{base}/starts"][s:e],
            "ends":    self.db[f"{base}/ends"][s:e],
            "widths":  self.db[f"{base}/widths"][s:e],
            "linkers": self.db[f"{base}/linkers"][s:e],
        }

    def get_methylation(
        self,
        fiber_id: str,
        chrom: str,
        mod_type: str = "5mC",
        sample: str | None = None,
    ) -> dict[str, np.ndarray]:
        """Return one fiber's modification calls of a given type.

        Works for 5mC, 5hmC, and 6mA. The returned dict always contains 'positions';
        for 5mC/5hmC it also contains 'probabilities' and 'is_methylated' if those
        are stored. 6mA has no per-call probability in our schema.

        Args:
            fiber_id: UUID string.
            chrom: Chromosome name.
            mod_type: One of "5mC", "5hmC", "6mA".
            sample: Sample name; defaults to the first sample.

        Returns:
            Dict of numpy arrays. Empty if the layer is missing or the fiber has
            no calls of this type.
        """
        sample = sample or self.samples[0]
        sl = self._get_fiber_slice(fiber_id, chrom, sample, mod_type)
        if sl is None:
            return {}
        s, e = sl
        base = schema.layer_path(sample, chrom, mod_type)
        result: dict[str, np.ndarray] = {
            "positions": self.db[f"{base}/positions"][s:e],
        }
        if f"{base}/probabilities" in self.db:
            result["probabilities"] = self.db[f"{base}/probabilities"][s:e]
        if f"{base}/is_methylated" in self.db:
            result["is_methylated"] = self.db[f"{base}/is_methylated"][s:e]
        return result

    def get_msp(
        self, fiber_id: str, chrom: str, sample: str | None = None
    ) -> dict[str, np.ndarray]:
        """Return one fiber's MSP (methylation-sensitive patch) intervals.

        Returns:
            Dict with keys 'starts', 'ends', 'widths'. Empty if the fiber has no
            MSPs or MSPs aren't stored for this chromosome.
        """
        sample = sample or self.samples[0]
        sl = self._get_fiber_slice(fiber_id, chrom, sample, "msp")
        if sl is None:
            return {}
        s, e = sl
        base = schema.layer_path(sample, chrom, "msp")
        return {
            "starts": self.db[f"{base}/starts"][s:e],
            "ends":   self.db[f"{base}/ends"][s:e],
            "widths": self.db[f"{base}/widths"][s:e],
        }

    def _find_fiber_chromosome(
        self, fiber_id: str, sample: str | None = None
    ) -> str | None:
        """Look up which chromosome a fiber lives on (from the per-sample fiber_lookup).

        Cached lazily. Returns None if the fiber isn't found in this sample.
        """
        sample = sample or self.samples[0]
        if sample not in self._fiber_lookup_cache:
            lookup = schema.fiber_lookup_path(sample)
            fids = self.db[f"{lookup}/fiber_ids"][:]
            chroms = self.db[f"{lookup}/chromosomes"][:]
            self._fiber_lookup_cache[sample] = dict(zip(fids, chroms, strict=False))
        r = self._fiber_lookup_cache[sample].get(fiber_id.encode())
        return r.decode() if r else None

    # ==================================================================
    # ANNOTATIONS + BULK QUERIES (Batch 3 - to come)
    # ==================================================================

    def query_annotation_fast(
        self,
        annotation_name: str,
        sample: str | None = None,
        feature_types: list[str] | None = None,
        max_regions: int | None = None,
    ) -> pd.DataFrame:
        """TODO: implement in Batch 3."""
        raise NotImplementedError("Will be added in Batch 3")
