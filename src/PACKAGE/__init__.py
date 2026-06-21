"""PACKAGE — single-molecule epigenomic analysis for long-read Fiber-seq data.

Public API:

    from PACKAGE.db import FiberDatabase
    from PACKAGE.db.builder import build_database
    from PACKAGE.extract import extract_ont
    from PACKAGE.viz import single_molecule_heatmap

See README.md for a quickstart, examples/ for tutorial notebooks, and docs/ for the
full API reference.
"""

__version__ = "0.4.0"

# Re-export the most commonly used names so users can do `from PACKAGE import X`.
# Kept minimal — most users should import from submodules.
from PACKAGE.db.database import FiberDatabase

__all__ = ["FiberDatabase", "__version__"]
