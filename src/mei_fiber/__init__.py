"""MEI-Fiber: Multimodal Epigenetic Integration with Fiber-seq.

Public API:

    from mei_fiber.db import FiberDatabase
    from mei_fiber.db.builder import build_database
    from mei_fiber.extract import extract_ont
    from mei_fiber.viz import single_molecule_heatmap

See README.md for a quickstart, examples/ for tutorial notebooks, and docs/ for the
full API reference.
"""

__version__ = "0.5.0"

# Re-export the most commonly used names so users can do `from mei_fiber import X`.
# Kept minimal — most users should import from submodules.
from mei_fiber.db.database import FiberDatabase

__all__ = ["FiberDatabase", "__version__"]
