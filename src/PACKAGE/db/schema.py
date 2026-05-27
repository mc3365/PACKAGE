"""Canonical HDF5 schema for PACKAGE databases.

This module is the **contract** between the ingestion (``PACKAGE.io``) and query
(``PACKAGE.db.database``) layers. Both sides import from here. To change the schema:
bump ``SCHEMA_VERSION``, update this module, and update both consumers in lockstep.

Layout
======

::

    /fiber_database.h5
    │
    ├── metadata                       (group with attrs: version, genome_version, ...)
    │
    ├── annotations/                   (shared across samples)
    │   ├── master_v3/
    │   │   └── features/<name>        (structured array: chr, start, end, id)
    │   └── FIRE_peaks/
    │       ├── <sample>_peaks
    │       └── <sample>_peak_annotations
    │
    └── <sample>/                      (one group per sample, e.g. 'd0', 'd4')
        ├── fiber_lookup/
        │   ├── fiber_ids              (sorted S50 array of UUID strings)
        │   └── chromosomes            (parallel S10 array)
        │
        └── <chrom>/                   (one per chromosome with data)
            ├── fiber_id_table         (S50 array; index = uint32 fiber ID)
            ├── fiber_metadata/
            │   ├── fiber_int_ids      (uint32)
            │   ├── starts             (uint32)
            │   └── ends               (uint32)
            ├── nucleosomes/
            │   ├── fiber_int_ids
            │   ├── starts
            │   ├── ends
            │   ├── widths
            │   └── linkers
            ├── 5mC/  6mA/  5hmC/      (each: fiber_int_ids, positions, probabilities, is_methylated)
            ├── msp/                   (fiber_int_ids, starts, ends, widths)
            └── _indices/              (slice-into-array lookup tables, keyed by uint32 fiber_int_id)
                ├── nucleosomes_slices
                ├── 5mC_slices
                └── ...

Design notes
============

- All coordinates are 0-based (consistent with BED).
- Within each chromosome, data arrays are sorted by ``fiber_int_id`` ascending. The
  ``_indices`` group gives O(1) row ranges per fiber.
- Per-chromosome integer IDs (not globally unique) — keeps IDs as ``uint32`` even for
  the largest datasets. The ``fiber_id_table`` translates back to the global UUID.
"""

from __future__ import annotations

import numpy as np

# ---------------------------------------------------------------------------
# Schema version
# ---------------------------------------------------------------------------

SCHEMA_VERSION: str = "8.0"
"""Bump on incompatible schema changes. Builders write this; readers check it."""


# ---------------------------------------------------------------------------
# Reserved top-level group names (not samples)
# ---------------------------------------------------------------------------

RESERVED_TOP_LEVEL: frozenset[str] = frozenset({"annotations", "metadata", "fiber_lookup"})
"""Names at the root that are NOT samples — used to filter when enumerating samples."""


# ---------------------------------------------------------------------------
# Modification types tracked per chromosome
# ---------------------------------------------------------------------------

MOD_TYPES: tuple[str, ...] = ("nucleosomes", "5mC", "5hmC", "6mA", "msp")
"""All per-fiber modification/feature groups in the chromosome layout."""


# ---------------------------------------------------------------------------
# Compound dtypes
# ---------------------------------------------------------------------------

SLICE_INDEX_DTYPE = np.dtype([
    ("fiber_int_id", np.uint32),
    ("start", np.uint32),
    ("end", np.uint32),
])
"""dtype for the per-fiber slice-into-array lookups stored in ``_indices/<mod>_slices``."""


ANNOTATION_REGION_DTYPE = np.dtype([
    ("chr", "S10"),
    ("start", np.uint32),
    ("end", np.uint32),
    ("id", "S50"),
])
"""dtype for annotation feature tables stored under ``annotations/master_v3/features/<name>``."""


FIRE_PEAK_DTYPE = np.dtype([
    ("fire_id", "S20"),
    ("chr", "S10"),
    ("start", np.uint32),
    ("end", np.uint32),
])
"""dtype for de-duplicated FIRE peak tables."""


FIRE_PEAK_ANNOTATION_DTYPE = np.dtype([
    ("fire_id", "S20"),
    ("annotation_id", "S50"),
])
"""dtype for FIRE-peak ↔ annotation many-to-many mapping table."""


# ---------------------------------------------------------------------------
# Path helpers — single source of truth for HDF5 paths
# ---------------------------------------------------------------------------

def fiber_metadata_path(sample: str, chrom: str) -> str:
    return f"{sample}/{chrom}/fiber_metadata"


def fiber_id_table_path(sample: str, chrom: str) -> str:
    return f"{sample}/{chrom}/fiber_id_table"


def mod_data_path(sample: str, chrom: str, mod: str) -> str:
    return f"{sample}/{chrom}/{mod}"


def slice_index_path(sample: str, chrom: str, mod: str) -> str:
    return f"{sample}/{chrom}/_indices/{mod}_slices"


def fiber_lookup_path(sample: str) -> str:
    return f"{sample}/fiber_lookup"


def annotation_features_path(annotation_set: str = "master_v3") -> str:
    return f"annotations/{annotation_set}/features"
