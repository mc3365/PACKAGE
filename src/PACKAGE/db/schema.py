"""Canonical HDF5 schema for PACKAGE databases.

This module is the **contract** between the build layer (``PACKAGE.db.builder``) and the
query layer (``PACKAGE.db.database``). Both sides import from here. To change the
schema: bump :data:`SCHEMA_VERSION`, update this module, and update both consumers in
lockstep.

The schema is a port of the V8 layout from ``build_fiber_database_v8.py`` /
``fiber_database_v8.py`` with three changes:

  1. All path literals are exposed as constants or path-builder functions instead of
     being hardcoded inline.
  2. All structured dtypes are defined once and imported everywhere they're used.
  3. The set of stored modification layers is **not fixed**: builders write whichever
     layers were extracted for a sample, and readers gracefully report layers as
     missing if not present. See :data:`SUPPORTED_LAYERS`.

Layout
======

::

    /database.h5
    |
    +-- metadata                       (group with attrs: version, genome_version, ...)
    |
    +-- annotations/                   (shared across samples)
    |   +-- master/
    |       +-- features/<type>        (structured array: chr, start, end, id)
    |                                  e.g. features/Compartment_B, features/Typical_Enhancer
    |
    +-- <sample>/                      (one group per sample; sample name is arbitrary)
        +-- fiber_lookup/
        |   +-- fiber_ids              (sorted S50 array of UUID strings)
        |   +-- chromosomes            (parallel S10 array)
        |
        +-- <chrom>/                   (one per chromosome with data)
            +-- fiber_id_table         (S50 array; row index = uint32 fiber ID)
            +-- fiber_metadata/
            |   +-- fiber_int_ids      (uint32)
            |   +-- starts             (uint32)
            |   +-- ends               (uint32)
            |
            +-- nucleosomes/           (present only if nucleosomes layer was extracted)
            |   +-- fiber_int_ids
            |   +-- starts
            |   +-- ends
            |   +-- widths
            |   +-- linkers
            |
            +-- 5mC/                   (each layer present only if extracted)
            |   +-- fiber_int_ids
            |   +-- positions
            |   +-- probabilities
            |   +-- is_methylated
            |
            +-- 5hmC/  6mA/            (same shape as 5mC; 6mA has no probabilities)
            |
            +-- msp/                   (fiber_int_ids, starts, ends, widths)
            |
            +-- fire_accessibility/    (fiber_int_ids, starts, ends, widths,
            |                           scores, haplotypes)
            |
            +-- _indices/              (slice-into-array lookup tables, keyed by uint32)
                +-- nucleosomes_slices
                +-- 5mC_slices
                +-- ...                (one per stored layer)

Design notes
============

- **Coordinates are 0-based** throughout (consistent with BED and modkit ref_position).
- **Per-chromosome integer fiber IDs** (not globally unique) keep IDs as ``uint32``
  even for the largest datasets. The ``fiber_id_table`` array translates back to the
  global UUID: index N gives the UUID for fiber_int_id N.
- **Data arrays are sorted by fiber_int_id ascending** within each chromosome. The
  ``_indices`` group gives O(1) row ranges per fiber via a structured array indexed
  by integer fiber ID.
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

RESERVED_TOP_LEVEL: frozenset[str] = frozenset({"annotations", "metadata"})
"""Names at the root that are NOT samples. Used to filter when enumerating samples."""


# ---------------------------------------------------------------------------
# Supported modification layers
# ---------------------------------------------------------------------------

SUPPORTED_LAYERS: tuple[str, ...] = (
    "5mC",
    "5hmC",
    "6mA",
    "msp",
    "nucleosomes",
    "fire_accessibility",
)
"""All modification/feature layers the schema knows about.

Each sample / chromosome stores whichever subset was extracted. Layer names here are
also the HDF5 group names under each chromosome.
"""

# Which arrays live inside each layer group. Used by builder for writing, by reader
# for validation. Keep in sync with the layout in this module's docstring.
LAYER_FIELDS: dict[str, tuple[str, ...]] = {
    "nucleosomes": ("fiber_int_ids", "starts", "ends", "widths", "linkers"),
    "5mC":         ("fiber_int_ids", "positions", "probabilities", "is_methylated"),
    "5hmC":        ("fiber_int_ids", "positions", "probabilities", "is_methylated"),
    "6mA":         ("fiber_int_ids", "positions"),
    "msp":         ("fiber_int_ids", "starts", "ends", "widths"),
    "fire_accessibility": (
        "fiber_int_ids",
        "starts",
        "ends",
        "widths",
        "scores",
        "haplotypes",
    ),
}


# ---------------------------------------------------------------------------
# Structured dtypes
# ---------------------------------------------------------------------------

SLICE_INDEX_DTYPE = np.dtype([
    ("fiber_int_id", np.uint32),
    ("start", np.uint32),
    ("end", np.uint32),
])
"""Per-fiber slice-into-array lookup stored in ``_indices/<layer>_slices``.

For each fiber, gives the [start, end) row range in the layer's data arrays. This is
how we go from "give me fiber X's 5mC calls" to a contiguous array slice in O(1).
"""


ANNOTATION_REGION_DTYPE = np.dtype([
    ("chr", "S10"),
    ("start", np.uint32),
    ("end", np.uint32),
    ("id", "S50"),
])
"""Annotation feature tables stored under ``annotations/master/features/<type>``.

The ``id`` field is the full annotation ID from the user's master BED file
(e.g. b"Typical_Enhancer_42"). The ``<type>`` group name is the ID prefix
(e.g. "Typical_Enhancer").
"""


# ---------------------------------------------------------------------------
# Path helpers - single source of truth for HDF5 paths
# ---------------------------------------------------------------------------

def metadata_path() -> str:
    """Path to the top-level metadata group."""
    return "metadata"


def fiber_lookup_path(sample: str) -> str:
    """Path to a sample's global fiber lookup group."""
    return f"{sample}/fiber_lookup"


def chrom_group_path(sample: str, chrom: str) -> str:
    """Path to a (sample, chromosome) group."""
    return f"{sample}/{chrom}"


def fiber_id_table_path(sample: str, chrom: str) -> str:
    """Path to the int->UUID lookup array for one (sample, chromosome)."""
    return f"{sample}/{chrom}/fiber_id_table"


def fiber_metadata_path(sample: str, chrom: str) -> str:
    """Path to the per-fiber metadata (starts, ends) group."""
    return f"{sample}/{chrom}/fiber_metadata"


def layer_path(sample: str, chrom: str, layer: str) -> str:
    """Path to a modification/feature layer group.

    Raises:
        ValueError: if ``layer`` is not in :data:`SUPPORTED_LAYERS`.
    """
    if layer not in SUPPORTED_LAYERS:
        raise ValueError(
            f"Unknown layer {layer!r}. Supported: {SUPPORTED_LAYERS}"
        )
    return f"{sample}/{chrom}/{layer}"


def slice_index_path(sample: str, chrom: str, layer: str) -> str:
    """Path to the per-fiber slice-into-array lookup for one layer."""
    if layer not in SUPPORTED_LAYERS:
        raise ValueError(
            f"Unknown layer {layer!r}. Supported: {SUPPORTED_LAYERS}"
        )
    return f"{sample}/{chrom}/_indices/{layer}_slices"


def annotation_features_path(annotation_type: str | None = None) -> str:
    """Path to the annotations group, or to a specific feature type within it.

    Args:
        annotation_type: e.g. "Typical_Enhancer". If None, returns the parent
            features/ group path.
    """
    if annotation_type is None:
        return "annotations/master/features"
    return f"annotations/master/features/{annotation_type}"


# ---------------------------------------------------------------------------
# Metadata attribute keys (write once at build, read by everyone)
# ---------------------------------------------------------------------------

# Stored under /metadata as HDF5 attrs. Defined here so builder and reader can't
# disagree on names.
META_SCHEMA_VERSION = "schema_version"
META_GENOME_VERSION = "genome_version"
META_METHYLATION_THRESHOLD = "methylation_threshold"
META_MIN_MSP_SIZE = "min_msp_size"
META_BUILD_DATE = "build_date"
META_PACKAGE_VERSION = "package_version"
