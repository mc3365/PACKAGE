"""Reproducibility manifest helpers for extraction runs."""

from __future__ import annotations

import hashlib
import json
import socket
import subprocess
from datetime import datetime
from pathlib import Path
from typing import Any

from PACKAGE import __version__


def write_manifest(
    output_dir: Path,
    manifest_data: dict[str, Any],
    filename: str = "PACKAGE_manifest.json",
) -> Path:
    """Write a human-readable JSON extraction manifest."""
    output_dir.mkdir(parents=True, exist_ok=True)
    payload = {
        "package_version": __version__,
        "created_at": datetime.now().isoformat(),
        "host": socket.gethostname(),
        **manifest_data,
    }
    path = output_dir / filename
    path.write_text(json.dumps(payload, indent=2, sort_keys=True, default=str) + "\n")
    return path


def capture_tool_versions(executables: dict[str, str]) -> dict[str, str]:
    """Capture the first nonempty output line from each tool's ``--version``."""
    versions: dict[str, str] = {}
    for name, executable in executables.items():
        result = subprocess.run(
            [executable, "--version"],
            check=False,
            capture_output=True,
            text=True,
        )
        combined = "\n".join((result.stdout, result.stderr))
        first_line = next((line.strip() for line in combined.splitlines() if line.strip()), "")
        versions[name] = first_line or f"unknown (exit {result.returncode})"
    return versions


def sha256_file(path: Path, chunk_size: int = 1 << 20) -> str:
    """Compute SHA256 without loading the whole file into memory."""
    digest = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()
