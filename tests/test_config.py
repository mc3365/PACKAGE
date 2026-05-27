"""Tests for PACKAGE.config."""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml

from PACKAGE.config import Config, load_config


def _write_yaml(path: Path, data: dict) -> Path:
    path.write_text(yaml.safe_dump(data))
    return path


def test_load_valid_config(tmp_path):
    """A well-formed YAML should parse into a Config instance."""
    cfg_path = _write_yaml(
        tmp_path / "config.yaml",
        {
            "output_dir": str(tmp_path / "out"),
            "output_file": "fibers.h5",
            "samples": [
                {
                    "name": "d0",
                    "input_dir": str(tmp_path),
                    "ft_extract": "d0.bed.gz",
                    "modkit_extract": "d0.mods.tsv.gz",
                }
            ],
            "annotations": {},
            "parameters": {"methylation_threshold": 0.5, "min_msp_size": 10},
        },
    )
    cfg = load_config(cfg_path)
    assert isinstance(cfg, Config)
    assert cfg.sample_names == ["d0"]
    assert cfg.parameters.methylation_threshold == 0.5
    assert cfg.output_path == tmp_path / "out" / "fibers.h5"


def test_load_missing_file_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        load_config(tmp_path / "no_such_file.yaml")


def test_invalid_threshold_rejected(tmp_path):
    """Pydantic should reject methylation_threshold > 1.0 at load time."""
    from pydantic import ValidationError

    cfg_path = _write_yaml(
        tmp_path / "bad.yaml",
        {
            "output_dir": str(tmp_path),
            "samples": [
                {"name": "d0", "input_dir": str(tmp_path), "ft_extract": "x.bed"}
            ],
            "parameters": {"methylation_threshold": 1.5},
        },
    )
    with pytest.raises(ValidationError):
        load_config(cfg_path)
