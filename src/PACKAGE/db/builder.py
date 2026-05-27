"""Database builder: extracted intermediate files → HDF5 database.

Entry point called by ``PACKAGE build``. Reads the intermediate files produced by
``PACKAGE extract`` (or by manual modkit/fibertools-rs invocations) and writes them
into the HDF5 schema defined in ``PACKAGE.db.schema``.

Responsibilities:
  1. Validate intermediate files exist
  2. Build per-chromosome string→int fiber ID maps
  3. Write data arrays + slice indices to HDF5 following schema
  4. Load shared annotations
  5. Write metadata attrs (version, parameters, build date)

STATUS: Stub. TODO(refactor): port from existing build_fiber_database_v8.py
"""

from __future__ import annotations

from PACKAGE.config import Config
from PACKAGE.utils import get_logger

log = get_logger(__name__)


def build_database(config: Config, samples: list[str] | None = None) -> None:
    """Build an HDF5 database from extracted intermediate files.

    Args:
        config: Validated PACKAGE config.
        samples: Subset of sample names. If None, builds all in config.
    """
    samples = samples or config.sample_names
    log.info(f"Building {config.output_path} for samples={samples}")
    # TODO(refactor): port logic from existing build_fiber_database_v8.py:
    #   - _validate_files
    #   - per-sample: read intermediates → assemble arrays → write HDF5
    #   - _load_annotations
    #   - write metadata
    raise NotImplementedError(
        "Builder orchestration to be ported from existing build_fiber_database_v8.py"
    )
