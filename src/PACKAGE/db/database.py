"""Read-only query interface for PACKAGE HDF5 databases.

This is the **library face** of the storage layer. It is a refactor of the existing
``fiber_database_v8.py`` from your scratch directory:

  - same public methods (``get_fibers_at``, ``get_methylation``, ``get_nucleosomes``,
    ``query_annotation_fast``, ``compare_samples``, ...)
  - imports schema paths from ``PACKAGE.db.schema`` instead of duplicating literals
  - spatial-index code split out into ``PACKAGE.db.spatial_index``

STATUS: Stub. Methods raise NotImplementedError until refactor lands.
TODO(refactor): port methods from uploads/fiber_database_v8.py here.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import h5py

from PACKAGE.db import schema
from PACKAGE.utils import get_logger

log = get_logger(__name__)


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

        # Cached lookups.
        self._chromosomes: dict[str, list[str]] = {}
        self._spatial_index: dict[str, dict[str, Any]] = {}
        self._int_map_cache: dict[tuple[str, str], dict[bytes, int]] = {}
        self._str_map_cache: dict[tuple[str, str], Any] = {}

        # Read metadata if present.
        if "metadata" in self.db:
            meta = self.db["metadata"]
            self.version: str = meta.attrs.get("version", "unknown")
            self.genome_version: str = meta.attrs.get("genome_version", "unknown")
            self.methylation_threshold: float = float(
                meta.attrs.get("methylation_threshold", 0.5)
            )
        else:
            self.version = "unknown"
            self.genome_version = "unknown"
            self.methylation_threshold = 0.5

        # TODO(refactor): port spatial index logic.
        if build_index:
            log.info("Spatial index build not yet ported — coming in next refactor pass.")

        log.info(f"FiberDatabase opened: {self.db_path.name}, samples={self.samples}")

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
    # Summary
    # ------------------------------------------------------------------
    def get_summary(self) -> dict[str, Any]:
        """Return a small dict describing the database contents."""
        return {
            "path": str(self.db_path),
            "version": self.version,
            "genome": self.genome_version,
            "samples": self.samples,
            "n_annotations": len(self.list_annotations()),
        }

    def get_chromosomes(self, sample: str | None = None) -> list[str]:
        """List chromosome names for a sample (filtered to those starting with 'chr')."""
        sample = sample or self.samples[0]
        if sample not in self._chromosomes:
            self._chromosomes[sample] = sorted(
                k for k in self.db[sample].keys() if k.startswith("chr")
            )
        return self._chromosomes[sample]

    def list_annotations(self) -> list[str]:
        """List available annotation sets under ``annotations/master_v3/features/``."""
        path = schema.annotation_features_path()
        if path in self.db:
            return sorted(self.db[path].keys())
        return []

    # ------------------------------------------------------------------
    # Query methods — TO BE PORTED FROM uploads/fiber_database_v8.py
    # ------------------------------------------------------------------
    def get_fibers_at(
        self, chrom: str, start: int, end: int, sample: str | None = None
    ) -> list[str]:
        """Return fiber IDs (as strings) overlapping ``[start, end)`` on ``chrom``."""
        raise NotImplementedError("Port from fiber_database_v8.FiberDatabase.get_fibers_at")

    def get_nucleosomes(
        self, fiber_id: str, chrom: str, sample: str | None = None
    ) -> dict[str, Any]:
        raise NotImplementedError("Port from fiber_database_v8.FiberDatabase.get_nucleosomes")

    def get_methylation(
        self,
        fiber_id: str,
        chrom: str,
        mod_type: str = "5mC",
        sample: str | None = None,
    ) -> dict[str, Any]:
        raise NotImplementedError("Port from fiber_database_v8.FiberDatabase.get_methylation")

    def get_msp(
        self, fiber_id: str, chrom: str, sample: str | None = None
    ) -> dict[str, Any]:
        raise NotImplementedError("Port from fiber_database_v8.FiberDatabase.get_msp")

    def query_annotation_fast(
        self,
        annotation_name: str,
        sample: str | None = None,
        feature_types: list[str] | None = None,
        max_regions: int | None = None,
    ) -> Any:  # pandas.DataFrame, untyped to avoid import at top
        raise NotImplementedError(
            "Port from fiber_database_v8.FiberDatabase.query_annotation_fast"
        )
