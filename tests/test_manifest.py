from __future__ import annotations

import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

from PACKAGE.extract.manifest import (
    capture_tool_versions,
    sha256_file,
    write_manifest,
)


def test_sha256_file_streams_expected_digest(tmp_path):
    path = tmp_path / "input.bin"
    path.write_bytes(b"fiber-seq\n" * 100)
    assert sha256_file(path, chunk_size=7) == hashlib.sha256(path.read_bytes()).hexdigest()


def test_write_manifest_adds_provenance_fields(tmp_path):
    path = write_manifest(
        tmp_path,
        {"sample": "d0", "input": Path("/data/d0.bam")},
        filename="PACKAGE_manifest_d0.json",
    )
    data = json.loads(path.read_text())
    assert data["sample"] == "d0"
    assert data["input"] == "/data/d0.bam"
    assert data["package_version"]
    assert data["created_at"]
    assert data["host"]


def test_capture_tool_versions_uses_first_nonempty_line(monkeypatch):
    def fake_run(command, **kwargs):
        assert command[1] == "--version"
        return SimpleNamespace(stdout="\nmodkit 0.5.0\nextra", stderr="", returncode=0)

    monkeypatch.setattr("PACKAGE.extract.manifest.subprocess.run", fake_run)
    versions = capture_tool_versions({"modkit": "/tools/modkit"})
    assert versions == {"modkit": "modkit 0.5.0"}
