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

    # ------------------------------------------------------------------
    # Class-level constants
    # ------------------------------------------------------------------
    MAX_FIBER_SPAN: int = 100_000
    """Fibers with ``end - start > MAX_FIBER_SPAN`` are treated as artifacts and
    excluded from all overlap queries. 100 kb is the V8 default; it filters out
    rare misaligned reads whose reported span is much larger than any plausible
    biological fiber. Filtering happens at the query layer, not the storage
    layer — the underlying HDF5 retains all fibers.
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
                    self.db, sample, self._get_id_maps,
                    max_fiber_span=self.MAX_FIBER_SPAN,
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

        Fibers whose reported span exceeds :attr:`MAX_FIBER_SPAN` are excluded
        from results (likely alignment artifacts).

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
        # The index was built with monster fibers already excluded (see spatial_index.py),
        # so we don't need to filter again here.
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
        # Filter monster fibers (artifact entries with implausibly large spans).
        spans = ends - starts
        valid = spans <= self.MAX_FIBER_SPAN
        starts, ends, fids = starts[valid], ends[valid], fids[valid]
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
        # Slice indices are written in ascending fiber_int_id order. Binary-search
        # the on-disk dataset so repeated per-fiber access does not load a complete
        # chromosome index for every lookup.
        idx = self.db[slice_path]
        target = int(int_id)
        low, high = 0, len(idx)
        while low < high:
            mid = (low + high) // 2
            row = idx[mid]
            row_id = int(row["fiber_int_id"])
            if row_id < target:
                low = mid + 1
            else:
                high = mid
        if low >= len(idx):
            return None
        row = idx[low]
        if int(row["fiber_int_id"]) != target:
            return None
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
    # BULK CHROMOSOME ACCESS (Batch 3 — load once, reuse for many regions)
    # ==================================================================
    # These methods support a workflow where a user wants to write their own
    # per-region loop without going through query_annotation_fast. Call
    # load_chromosome_data() once per chromosome, then use the returned dict
    # with get_overlapping_fiber_ids / get_spanning_fiber_ids for fast,
    # numpy-only overlap selection.

    def load_chromosome_data(
        self,
        sample: str,
        chrom: str,
        layers: list[str] | None = None,
    ) -> dict[str, Any] | None:
        """Load all arrays for one chromosome into memory.

        Call once per chromosome, then reuse for many region queries via
        :meth:`get_overlapping_fiber_ids` / :meth:`get_spanning_fiber_ids`.
        Monster fibers (span > :attr:`MAX_FIBER_SPAN`) are excluded automatically
        from the fiber metadata arrays — but NOT from per-layer arrays, so
        per-fiber lookups by integer ID still work.

        Args:
            sample: Sample name.
            chrom: Chromosome name.
            layers: Layers to load. Default: ``["nucleosomes"]``. Pass
                ``["nucleosomes", "5mC", "5hmC", "6mA", "msp"]`` to load all.

        Returns:
            A dict with these top-level keys:

            - ``fm_ids``, ``fm_starts``, ``fm_ends``: filtered fiber metadata arrays
            - ``id_table``: the raw fiber_id_table (uint32 -> UUID bytestring)
            - ``<layer>``: ``{field: ndarray}`` for each loaded layer
            - ``<layer>_idx``: ``{fiber_int_id: (start_row, end_row)}`` slice lookup

            Returns ``None`` if the chromosome isn't in the database.
        """
        base = schema.chrom_group_path(sample, chrom)
        if base not in self.db:
            return None

        cd: dict[str, Any] = {}

        # Fiber metadata (will be filtered).
        meta = schema.fiber_metadata_path(sample, chrom)
        cd["fm_ids"] = self.db[f"{meta}/fiber_int_ids"][:]
        cd["fm_starts"] = self.db[f"{meta}/starts"][:]
        cd["fm_ends"] = self.db[f"{meta}/ends"][:]

        # Filter monster fibers and log how many were dropped.
        spans = cd["fm_ends"] - cd["fm_starts"]
        valid = spans <= self.MAX_FIBER_SPAN
        n_excluded = int((~valid).sum())
        if n_excluded > 0:
            n_before = len(cd["fm_ids"])
            cd["fm_ids"] = cd["fm_ids"][valid]
            cd["fm_starts"] = cd["fm_starts"][valid]
            cd["fm_ends"] = cd["fm_ends"][valid]
            log.info(
                f"  {chrom}: excluded {n_excluded:,}/{n_before:,} monster fibers "
                f"(>{self.MAX_FIBER_SPAN/1000:.0f} kb)"
            )

        # ID lookup table (unfiltered — needed to translate any int ID to UUID).
        cd["id_table"] = self.db[schema.fiber_id_table_path(sample, chrom)][:]

        layers = list(layers) if layers else ["nucleosomes"]
        for layer in layers:
            layer_grp = schema.layer_path(sample, chrom, layer)
            if layer_grp not in self.db:
                continue
            lg = self.db[layer_grp]
            ld: dict[str, Any] = {}
            ld["fiber_int_ids"] = lg["fiber_int_ids"][:]

            if layer == "nucleosomes":
                ld["starts"] = lg["starts"][:]
                ld["ends"] = lg["ends"][:]
                ld["widths"] = lg["widths"][:]
                ld["linkers"] = lg["linkers"][:]
            elif layer in ("5mC", "5hmC"):
                ld["positions"] = lg["positions"][:]
                if "is_methylated" in lg:
                    ld["is_methylated"] = lg["is_methylated"][:]
                if "probabilities" in lg:
                    ld["probabilities"] = lg["probabilities"][:]
            elif layer == "6mA":
                ld["positions"] = lg["positions"][:]
            elif layer == "msp":
                ld["starts"] = lg["starts"][:]
                ld["ends"] = lg["ends"][:]
                ld["widths"] = lg["widths"][:]

            cd[layer] = ld

            # Slice index: fiber_int_id -> (start_idx, end_idx).
            idx_path = schema.slice_index_path(sample, chrom, layer)
            if idx_path in self.db:
                idx = self.db[idx_path][:]
                cd[f"{layer}_idx"] = {
                    int(row["fiber_int_id"]): (int(row["start"]), int(row["end"]))
                    for row in idx
                }
            else:
                cd[f"{layer}_idx"] = {}

        return cd

    @staticmethod
    def get_overlapping_fiber_ids(
        chrom_data: dict[str, Any], start: int, end: int
    ) -> np.ndarray:
        """Return fiber int IDs whose span overlaps ``[start, end)`` (V8 semantics).

        Pure numpy — no HDF5 I/O. Use after :meth:`load_chromosome_data`.

        Uses V8's inclusive overlap (``fm_ends >= start`` AND ``fm_starts <= end``)
        for output parity with :meth:`query_annotation_fast`.
        """
        mask = (
            (chrom_data["fm_ends"] >= start)
            & (chrom_data["fm_starts"] <= end)
        )
        return chrom_data["fm_ids"][mask]

    @staticmethod
    def get_spanning_fiber_ids(
        chrom_data: dict[str, Any],
        start: int,
        end: int,
        min_coverage: float = 0.8,
    ) -> np.ndarray:
        """Return fiber int IDs that cover at least ``min_coverage`` of ``[start, end)``.

        Pure numpy — no HDF5 I/O. Use after :meth:`load_chromosome_data`.
        Useful for analyses that need fibers spanning most of a regulatory region
        rather than barely touching it.

        Args:
            chrom_data: Output of :meth:`load_chromosome_data`.
            start, end: Region bounds (0-based, half-open).
            min_coverage: Minimum fraction of region length the fiber must cover.

        Returns:
            uint32 array of fiber int IDs.
        """
        reg_len = end - start
        if reg_len <= 0:
            return np.array([], dtype=np.uint32)
        mask = (
            (chrom_data["fm_ends"] >= start)
            & (chrom_data["fm_starts"] <= end)
        )
        if not mask.any():
            return np.array([], dtype=np.uint32)
        fids = chrom_data["fm_ids"][mask]
        starts = chrom_data["fm_starts"][mask]
        ends = chrom_data["fm_ends"][mask]
        overlap_start = np.maximum(starts, start)
        overlap_end = np.minimum(ends, end)
        coverage = np.maximum(overlap_end - overlap_start, 0) / reg_len
        return fids[coverage >= min_coverage]

    # ==================================================================
    # ANNOTATIONS (Batch 3)
    # ==================================================================

    def get_annotation_regions(
        self, annotation_name: str
    ) -> list[tuple[str, int, int]]:
        """Return all regions for a given annotation type as (chrom, start, end) tuples.

        Args:
            annotation_name: Annotation type prefix, e.g. "Typical_Enhancer", "CGI".

        Returns:
            List of (chrom, start, end) tuples. Empty if the annotation type doesn't
            exist in this database.
        """
        root = self._annotation_features_root()
        if root is None:
            return []
        path = f"{root}/{annotation_name}"
        if path not in self.db:
            available = self.list_annotations()
            log.info(
                f"Annotation type {annotation_name!r} not found. "
                f"Available (first 10): {available[:10]}"
            )
            return []
        data = self.db[path][:]
        return [
            (r["chr"].decode(), int(r["start"]), int(r["end"])) for r in data
        ]

    def get_annotation_by_id(
        self, annotation_id: str
    ) -> tuple[str, int, int] | None:
        """Look up a single annotation region by its full ID string.

        Annotation IDs encode the type as the prefix (e.g. "Typical_Enhancer_42"
        has type "Typical_Enhancer" and numeric suffix "42"). We strip the last
        underscore-separated token to get the type, then linear-scan that type's
        feature table for the matching ID.

        Args:
            annotation_id: Full annotation ID, e.g. "Typical_Enhancer_42".

        Returns:
            (chrom, start, end) tuple, or None if no match.
        """
        # Type = everything before the last underscore-separated token.
        anno_type = "_".join(annotation_id.split("_")[:-1])
        root = self._annotation_features_root()
        if root is None:
            return None
        path = f"{root}/{anno_type}"
        if path not in self.db:
            return None
        data = self.db[path][:]
        if "id" not in data.dtype.names:
            return None
        for row in data:
            if row["id"].decode() == annotation_id:
                return (row["chr"].decode(), int(row["start"]), int(row["end"]))
        return None

    def get_annotation_regions_with_ids(
        self, annotation_name: str
    ) -> list[tuple[str, int, int, str]]:
        """Same as :meth:`get_annotation_regions` but also returns each region's ID.

        Useful when downstream analyses need to refer back to specific annotation
        entries by their full ID (e.g. ``"Typical_Enhancer_42"``).

        Args:
            annotation_name: Annotation type prefix.

        Returns:
            List of ``(chrom, start, end, unique_id)`` tuples. Empty if the
            annotation type doesn't exist in this database.
        """
        root = self._annotation_features_root()
        if root is None:
            return []
        path = f"{root}/{annotation_name}"
        if path not in self.db:
            available = self.list_annotations()
            log.info(
                f"Annotation type {annotation_name!r} not found. "
                f"Available (first 10): {available[:10]}"
            )
            return []
        data = self.db[path][:]
        regions = []
        for row in data:
            if "id" in data.dtype.names:
                uid = row["id"].decode()
            else:
                uid = f"{annotation_name}_{len(regions)}"
            regions.append(
                (row["chr"].decode(), int(row["start"]), int(row["end"]), uid)
            )
        return regions

    @staticmethod
    def group_regions_by_chrom(
        regions: list[tuple[str, int, int, str]],
    ) -> dict[str, list[tuple[int, int, str]]]:
        """Group ``(chr, start, end, id)`` tuples by chromosome.

        Helper for users building their own per-chromosome loops on top of
        :meth:`load_chromosome_data`.
        """
        from collections import defaultdict

        by_chrom: dict[str, list[tuple[int, int, str]]] = defaultdict(list)
        for chrom, start, end, uid in regions:
            by_chrom[chrom].append((start, end, uid))
        return dict(by_chrom)

    # ==================================================================
    # BULK QUERIES (Batch 3)
    # ==================================================================

    # Default set of feature types if the caller doesn't specify. Matches V8.
    _DEFAULT_FEATURE_TYPES: tuple[str, ...] = ("nucleosomes", "5mC")

    def query_annotation_fast(
        self,
        annotation_name: str,
        sample: str | None = None,
        feature_types: list[str] | None = None,
        max_regions: int | None = None,
    ) -> pd.DataFrame:
        """Vectorized bulk query: per-fiber stats for each region of an annotation type.

        For each region of the given annotation type, finds all fibers overlapping it
        and computes per-fiber summary statistics for each requested layer:

        - **nucleosomes**: ``n_nucleosomes``, ``mean_nuc_width``, ``mean_linker``,
          ``mean_repeat_length``
        - **5mC**: ``n_cpg``, ``pct_methylated`` (or ``mean_methylation`` if
          ``is_methylated`` isn't stored)
        - **5hmC**: ``n_5hmC``, ``pct_hydroxymethylated``
        - **6mA**: ``n_6mA``
        - **msp**: ``n_msp``, ``mean_msp_width``

        The "fast" in the name refers to vectorization: each chromosome's data arrays
        are loaded once, then all regions on that chromosome are processed against
        them. This is dramatically faster than the per-fiber, per-region nested-loop
        approach for genome-scale queries.

        Position filters use V8's inclusive ``<=`` semantics (kept for output parity).

        Args:
            annotation_name: Annotation type prefix (e.g. "CGI", "Typical_Enhancer").
            sample: Sample to query. Defaults to the first sample.
            feature_types: Which layers to compute. Defaults to ["nucleosomes", "5mC"].
            max_regions: If set, only the first N regions are processed (useful for
                testing).

        Returns:
            DataFrame with one row per (region, fiber) combination. Empty DataFrame
            if the annotation type doesn't exist or there are no overlapping fibers.
        """
        sample = sample or self.samples[0]
        feature_types = list(feature_types) if feature_types else list(
            self._DEFAULT_FEATURE_TYPES
        )

        regions = self.get_annotation_regions(annotation_name)
        if not regions:
            return pd.DataFrame()
        if max_regions:
            regions = regions[:max_regions]
        log.info(
            f"Bulk query: {len(regions)} {annotation_name!r} regions, "
            f"sample={sample}, layers={feature_types}"
        )

        # Group regions by chromosome so we can load each chromosome's arrays once.
        from collections import defaultdict
        by_chrom: dict[str, list[tuple[int, int]]] = defaultdict(list)
        for c, s, e in regions:
            by_chrom[c].append((s, e))

        results: list[dict[str, Any]] = []
        for chrom, chrom_regions in by_chrom.items():
            base = schema.chrom_group_path(sample, chrom)
            if base not in self.db:
                continue
            _, i2s = self._get_id_maps(sample, chrom)

            # Load fiber metadata once for this chromosome.
            meta = schema.fiber_metadata_path(sample, chrom)
            fm_ids = self.db[f"{meta}/fiber_int_ids"][:]
            fm_starts = self.db[f"{meta}/starts"][:]
            fm_ends = self.db[f"{meta}/ends"][:]
            # Filter monster fibers (artifact entries with implausibly large spans).
            spans = fm_ends - fm_starts
            valid = spans <= self.MAX_FIBER_SPAN
            fm_ids, fm_starts, fm_ends = fm_ids[valid], fm_starts[valid], fm_ends[valid]

            # Load per-layer data arrays + slice lookups once per chromosome.
            layer_data, layer_slices = self._load_chrom_layer_data(
                sample, chrom, feature_types
            )

            # Process every region against the pre-loaded data.
            for reg_start, reg_end in chrom_regions:
                region_id = f"{chrom}:{reg_start}-{reg_end}"
                overlap_mask = (fm_ends >= reg_start) & (fm_starts <= reg_end)

                for fid_int in fm_ids[overlap_mask]:
                    rec: dict[str, Any] = {
                        "region_id": region_id,
                        "annotation": annotation_name,
                        "fiber_id": i2s[fid_int].decode(),
                        "sample": sample,
                    }
                    self._compute_per_fiber_stats(
                        rec, fid_int, reg_start, reg_end,
                        feature_types, layer_data, layer_slices,
                    )
                    results.append(rec)

        log.info(f"Extracted {len(results):,} (fiber, region) records")
        return pd.DataFrame(results)

    # alias preserved from V8
    def query_annotation(
        self,
        annotation_name: str,
        sample: str | None = None,
        feature_types: list[str] | None = None,
        max_regions: int | None = None,
    ) -> pd.DataFrame:
        """Alias of :meth:`query_annotation_fast` for V8 API compatibility."""
        return self.query_annotation_fast(
            annotation_name, sample, feature_types, max_regions
        )

    def _load_chrom_layer_data(
        self,
        sample: str,
        chrom: str,
        feature_types: list[str],
    ) -> tuple[dict[str, dict[str, np.ndarray]], dict[str, dict[int, tuple[int, int]]]]:
        """Load all layer arrays + slice-index dicts for one chromosome.

        Returns:
            (layer_data, layer_slices) where:

            - layer_data maps layer_name -> {field_name: numpy_array}, only including
              layers that are both requested AND present in the database.
            - layer_slices maps layer_name -> {fiber_int_id: (start_row, end_row)},
              for fast per-fiber row-range lookup.
        """
        base = schema.chrom_group_path(sample, chrom)
        layer_data: dict[str, dict[str, np.ndarray]] = {}

        # Each layer has a different set of fields; the schema's LAYER_FIELDS map
        # tells us what to load. We only load layers the caller asked for and that
        # are actually present in the database.
        for layer in feature_types:
            if layer not in schema.SUPPORTED_LAYERS:
                continue
            layer_grp = f"{base}/{layer}"
            if layer_grp not in self.db:
                continue
            fields_to_load = schema.LAYER_FIELDS[layer]
            layer_data[layer] = {
                field: self.db[f"{layer_grp}/{field}"][:]
                for field in fields_to_load
                if f"{layer_grp}/{field}" in self.db
            }

        # Load per-layer slice indices: layer -> {fiber_int_id: (start, end)}.
        # Doing this as a dict (rather than the structured array) makes per-fiber
        # lookups O(1) inside the region loop.
        layer_slices: dict[str, dict[int, tuple[int, int]]] = {}
        for layer in feature_types:
            if layer not in schema.SUPPORTED_LAYERS:
                continue
            slice_path = schema.slice_index_path(sample, chrom, layer)
            if slice_path not in self.db:
                continue
            idx = self.db[slice_path][:]
            layer_slices[layer] = {
                int(row["fiber_int_id"]): (int(row["start"]), int(row["end"]))
                for row in idx
            }

        return layer_data, layer_slices

    def _compute_per_fiber_stats(
        self,
        rec: dict[str, Any],
        fid_int: int,
        reg_start: int,
        reg_end: int,
        feature_types: list[str],
        layer_data: dict[str, dict[str, np.ndarray]],
        layer_slices: dict[str, dict[int, tuple[int, int]]],
    ) -> None:
        """Populate ``rec`` in-place with per-fiber stats for one (fiber, region).

        Each layer has its own field-name conventions in the output record, matching
        V8 exactly (n_nucleosomes, n_cpg, pct_methylated, etc.) so existing
        downstream code keeps working.

        Position-filter semantics use V8's inclusive ``<=`` (kept for output parity).
        """
        # ---- nucleosomes ----
        if "nucleosomes" in feature_types and "nucleosomes" in layer_data:
            sl = layer_slices.get("nucleosomes", {}).get(int(fid_int))
            if sl:
                s, e = sl
                d = layer_data["nucleosomes"]
                m = (d["starts"][s:e] >= reg_start) & (d["ends"][s:e] <= reg_end)
                rec["n_nucleosomes"] = int(m.sum())
                if m.any():
                    rec["mean_nuc_width"] = float(d["widths"][s:e][m].mean())
                    lk = d["linkers"][s:e][m]
                    lk = lk[lk >= 0]  # filter the -1 sentinel for missing linker
                    if len(lk) > 0:
                        rec["mean_linker"] = float(lk.mean())
                        rec["mean_repeat_length"] = (
                            rec["mean_nuc_width"] + rec["mean_linker"]
                        )
            else:
                rec["n_nucleosomes"] = 0

        # ---- 5mC ----
        if "5mC" in feature_types and "5mC" in layer_data:
            sl = layer_slices.get("5mC", {}).get(int(fid_int))
            if sl:
                s, e = sl
                d = layer_data["5mC"]
                m = (d["positions"][s:e] >= reg_start) & (d["positions"][s:e] <= reg_end)
                rec["n_cpg"] = int(m.sum())
                if m.any():
                    if "is_methylated" in d:
                        rec["pct_methylated"] = float(
                            d["is_methylated"][s:e][m].mean() * 100
                        )
                    elif "probabilities" in d:
                        rec["mean_methylation"] = float(
                            d["probabilities"][s:e][m].mean()
                        )
            else:
                rec["n_cpg"] = 0

        # ---- 5hmC ----
        if "5hmC" in feature_types and "5hmC" in layer_data:
            sl = layer_slices.get("5hmC", {}).get(int(fid_int))
            if sl:
                s, e = sl
                d = layer_data["5hmC"]
                m = (d["positions"][s:e] >= reg_start) & (d["positions"][s:e] <= reg_end)
                rec["n_5hmC"] = int(m.sum())
                if m.any() and "is_methylated" in d:
                    rec["pct_hydroxymethylated"] = float(
                        d["is_methylated"][s:e][m].mean() * 100
                    )
            else:
                rec["n_5hmC"] = 0

        # ---- 6mA ----
        if "6mA" in feature_types and "6mA" in layer_data:
            sl = layer_slices.get("6mA", {}).get(int(fid_int))
            if sl:
                s, e = sl
                d = layer_data["6mA"]
                m = (d["positions"][s:e] >= reg_start) & (d["positions"][s:e] <= reg_end)
                rec["n_6mA"] = int(m.sum())
            else:
                rec["n_6mA"] = 0

        # ---- msp ----
        if "msp" in feature_types and "msp" in layer_data:
            sl = layer_slices.get("msp", {}).get(int(fid_int))
            if sl:
                s, e = sl
                d = layer_data["msp"]
                m = (d["starts"][s:e] >= reg_start) & (d["ends"][s:e] <= reg_end)
                rec["n_msp"] = int(m.sum())
                if m.any():
                    rec["mean_msp_width"] = float(d["widths"][s:e][m].mean())
            else:
                rec["n_msp"] = 0

    def query_by_id(
        self,
        region_id: str,
        sample: str | None = None,
        feature_types: list[str] | None = None,
    ) -> pd.DataFrame:
        """Per-fiber query for a single annotation region, looked up by ID.

        Useful for inspecting one specific region without running a bulk query. Slower
        than :meth:`query_annotation_fast` per-region but trivial for a single lookup.

        Args:
            region_id: Annotation ID, e.g. "Typical_Enhancer_42".
            sample: Sample to query.
            feature_types: Layers to compute; defaults to all five supported layers.

        Returns:
            DataFrame with one row per fiber overlapping the region. Empty DataFrame
            if the region ID isn't found or no fibers overlap.
        """
        sample = sample or self.samples[0]
        region = self.get_annotation_by_id(region_id)
        if not region:
            log.info(f"Region {region_id!r} not found")
            return pd.DataFrame()
        chrom, start, end = region
        feature_types = list(feature_types) if feature_types else list(
            schema.SUPPORTED_LAYERS
        )

        fibers = self.get_fibers_at(chrom, start, end, sample)
        if not fibers:
            return pd.DataFrame()

        results: list[dict[str, Any]] = []
        for fid in fibers:
            rec: dict[str, Any] = {
                "region_id": region_id,
                "fiber_id": fid,
                "sample": sample,
            }

            if "nucleosomes" in feature_types:
                nucs = self.get_nucleosomes(fid, chrom, sample)
                if nucs:
                    m = (nucs["starts"] >= start) & (nucs["ends"] <= end)
                    rec["n_nucleosomes"] = int(m.sum())
                    if m.any():
                        rec["mean_nuc_width"] = float(nucs["widths"][m].mean())
                else:
                    rec["n_nucleosomes"] = 0

            if "5mC" in feature_types:
                mc = self.get_methylation(fid, chrom, "5mC", sample)
                if mc and len(mc["positions"]) > 0:
                    m = (mc["positions"] >= start) & (mc["positions"] <= end)
                    rec["n_cpg"] = int(m.sum())
                    if m.any() and "is_methylated" in mc:
                        rec["pct_methylated"] = float(
                            mc["is_methylated"][m].mean() * 100
                        )
                else:
                    rec["n_cpg"] = 0

            if "6mA" in feature_types:
                ma = self.get_methylation(fid, chrom, "6mA", sample)
                if ma and len(ma["positions"]) > 0:
                    # Note: V8 had a parenthesization bug here that would have
                    # raised at runtime. Fixed to int(boolean.sum()).
                    m = (ma["positions"] >= start) & (ma["positions"] <= end)
                    rec["n_6mA"] = int(m.sum())
                else:
                    rec["n_6mA"] = 0

            if "msp" in feature_types:
                msp = self.get_msp(fid, chrom, sample)
                if msp and len(msp["starts"]) > 0:
                    m = (msp["starts"] >= start) & (msp["ends"] <= end)
                    rec["n_msp"] = int(m.sum())
                else:
                    rec["n_msp"] = 0

            results.append(rec)

        return pd.DataFrame(results)

    # ==================================================================
    # CROSS-SAMPLE COMPARISON (Batch 3)
    # ==================================================================

    def compare_samples(
        self,
        annotation_name: str,
        feature_types: list[str] | None = None,
        max_regions: int | None = None,
    ) -> pd.DataFrame:
        """Run :meth:`query_annotation_fast` on every sample and concatenate.

        The returned DataFrame has the same per-(region, fiber) row structure as a
        single-sample query, with the additional ``sample`` column distinguishing
        which sample each row came from. This is the natural input to a comparison
        analysis (e.g. d0 vs d4, or ONT vs PacBio).

        Args:
            annotation_name: Annotation type prefix.
            feature_types: Layers to compute.
            max_regions: Optional per-sample region cap.

        Returns:
            Concatenated DataFrame across all samples in the database. Empty if no
            sample produced any results.
        """
        frames: list[pd.DataFrame] = []
        for sample in self.samples:
            log.info(f"--- {sample} ---")
            df = self.query_annotation_fast(
                annotation_name,
                sample=sample,
                feature_types=feature_types,
                max_regions=max_regions,
            )
            if not df.empty:
                frames.append(df)
        if not frames:
            return pd.DataFrame()
        combined = pd.concat(frames, ignore_index=True)
        log.info(
            f"Combined: {len(combined):,} records across {len(self.samples)} samples"
        )
        return combined
