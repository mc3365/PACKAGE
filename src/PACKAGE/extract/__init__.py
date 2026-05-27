"""Extraction layer: BAM → intermediate per-fiber tables.

Architecture
============

Platform-specific extraction is isolated in ``ont.py`` and ``pacbio.py``. Both modules
expose ``extract_samples(config, samples=...)`` with the same signature; the CLI
dispatches based on ``--platform``.

Each platform module orchestrates established tools:

- **ONT** (``ont.py``): modkit extract (5mC/5hmC) + fibertools-rs ft extract (m6A, MSP,
  nucleosomes)
- **PacBio** (``pacbio.py``): pb-CpG-tools or fibertools (5mC) + fibertools-rs ft extract
  (m6A, MSP, nucleosomes)

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
#     from PACKAGE.extract import extract_ont, extract_pacbio
# without remembering the submodule paths.
from PACKAGE.extract.ont import extract_samples as extract_ont
from PACKAGE.extract.pacbio import extract_samples as extract_pacbio

__all__ = ["extract_ont", "extract_pacbio"]
