"""HDF5 storage layer for PACKAGE.

Submodules:
    schema          Canonical HDF5 schema definitions (the contract between extract and db).
    database        FiberDatabase: read-only query interface.
    builder         build_database(): write extracted intermediates into HDF5.
    spatial_index   IntervalTree-based spatial index for fast region queries.
"""

from PACKAGE.db.database import FiberDatabase

__all__ = ["FiberDatabase"]
