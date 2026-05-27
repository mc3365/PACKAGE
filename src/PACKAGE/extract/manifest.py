"""Reproducibility manifest writer.

For every extraction run, write a JSON manifest capturing:

  - Tool versions (modkit, fibertools-rs, pb-CpG-tools, PACKAGE itself)
  - Parameters used
  - Input file paths and content hashes (e.g. SHA256 of BAM)
  - Output file paths
  - Run date and host

This manifest ships alongside the database in the supplementary materials and lets
reviewers / collaborators trace any analysis result back to its full provenance.

STATUS: Stub.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any


def write_manifest(output_dir: Path, manifest_data: dict[str, Any]) -> Path:
    """Write a JSON manifest to ``output_dir / 'PACKAGE_manifest.json'``."""
    raise NotImplementedError


def capture_tool_versions() -> dict[str, str]:
    """Run --version on each external tool and capture output.

    Returns:
        Dict of tool_name → version string.
    """
    raise NotImplementedError


def sha256_file(path: Path, chunk_size: int = 1 << 20) -> str:
    """Compute SHA256 of a file by streaming chunks."""
    raise NotImplementedError
