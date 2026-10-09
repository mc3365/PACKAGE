"""Extraction layer: BAM → intermediate per-fiber tables.

Architecture
============

Platform-specific extraction is isolated in ``ont.py`` and ``pacbio.py``. ONT is the
validated end-to-end implementation. PacBio Phase 1 normalizes fibertools outputs
into the same builder inputs used by the HDF5 schema.

Each platform module orchestrates established tools:

- **ONT** (``ont.py``): modkit extract (5mC/5hmC) + fibertools-rs ft extract (m6A, MSP,
  nucleosomes)
- **PacBio** (``pacbio.py``): fibertools-rs ``ft extract`` normalization for
  5mC, 6mA, MSP, and nucleosomes

Shared helpers live in:

- ``qc.py``: input validation and per-sample QC metrics
- ``manifest.py``: reproducibility manifest writer (tool versions, parameters, input
  hashes)

The extraction layer **wraps** these established tools rather than reimplementing them.
This is an intentional design choice — see manuscript Methods for rationale.

Submodules:
    ont          ONT extraction orchestration
    pacbio       PacBio extraction orchestration
    qc           Input validation + QC metrics
    manifest     Reproducibility manifest
"""

# Re-export the top-level entry points so users can do:
#     from mei_fiber.extract import extract_ont, extract_pacbio
# without remembering the submodule paths.
from mei_fiber.extract.ont import extract_samples as extract_ont
from mei_fiber.extract.pacbio import extract_samples as extract_pacbio

__all__ = ["extract_ont", "extract_pacbio"]
