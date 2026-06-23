"""Tests for PACKAGE.config."""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml

from PACKAGE.config import Config, load_config


def _write_yaml(path: Path, data: dict) -> Path:
    path.write_text(yaml.safe_dump(data))
    return path


def _minimal_sample_dict(tmp_path: Path, name: str = "d0") -> dict:
    """Return a minimal valid sample config dict."""
    return {
        "name": name,
        "layers": {
            "nucleosomes": str(tmp_path / f"{name}_nuc.csv"),
            "5mC":         str(tmp_path / f"{name}_raw_mods.tsv.gz"),
        },
    }


def test_load_valid_config(tmp_path):
    """A well-formed YAML should parse into a Config instance."""
    cfg_path = _write_yaml(
        tmp_path / "config.yaml",
        {
            "output_dir": str(tmp_path / "out"),
            "output_file": "fibers.h5",
            "samples": [_minimal_sample_dict(tmp_path, "d0")],
            "parameters": {"methylation_threshold": 0.5, "min_msp_size": 10},
        },
    )
    cfg = load_config(cfg_path)
    assert isinstance(cfg, Config)
    assert cfg.sample_names == ["d0"]
    assert cfg.parameters.methylation_threshold == 0.5
    assert cfg.extraction.threads == 16
    assert cfg.extraction.overwrite is False
    assert cfg.extraction.validate_bam is True
    assert cfg.extraction.write_manifest is True
    assert cfg.build.build_spatial_index is False
    assert cfg.output_path == tmp_path / "out" / "fibers.h5"


def test_build_settings_are_typed_and_validated(tmp_path):
    """Build settings should be recorded in the YAML schema."""
    cfg_path = _write_yaml(
        tmp_path / "build.yaml",
        {
            "output_dir": str(tmp_path),
            "samples": [_minimal_sample_dict(tmp_path)],
            "build": {"build_spatial_index": True},
        },
    )
    cfg = load_config(cfg_path)
    assert cfg.build.build_spatial_index is True


def test_extraction_settings_are_typed_and_validated(tmp_path):
    """Extraction settings should be configurable without affecting old configs."""
    cfg_path = _write_yaml(
        tmp_path / "extract.yaml",
        {
            "output_dir": str(tmp_path),
            "samples": [_minimal_sample_dict(tmp_path)],
            "extraction": {
                "threads": 8,
                "overwrite": True,
                "keep_nucleosome_bed": False,
                "validate_bam": False,
                "write_manifest": False,
            },
        },
    )
    cfg = load_config(cfg_path)
    assert cfg.extraction.threads == 8
    assert cfg.extraction.overwrite is True
    assert cfg.extraction.keep_nucleosome_bed is False
    assert cfg.extraction.validate_bam is False
    assert cfg.extraction.write_manifest is False


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
            "samples": [_minimal_sample_dict(tmp_path)],
            "parameters": {"methylation_threshold": 1.5},
        },
    )
    with pytest.raises(ValidationError):
        load_config(cfg_path)


# ---------------------------------------------------------------------------
# Per-layer file-path schema
# ---------------------------------------------------------------------------
def test_missing_nucleosomes_layer_rejected(tmp_path):
    """Nucleosomes is required; sample without it should fail validation."""
    from pydantic import ValidationError

    cfg_path = _write_yaml(
        tmp_path / "no_nuc.yaml",
        {
            "output_dir": str(tmp_path),
            "samples": [
                {
                    "name": "d0",
                    "layers": {
                        "5mC": str(tmp_path / "d0_raw_mods.tsv.gz"),
                    },
                },
            ],
        },
    )
    with pytest.raises(ValidationError, match="nucleosomes"):
        load_config(cfg_path)


def test_unknown_layer_name_rejected(tmp_path):
    """Unknown layer names should be caught at config-load time."""
    from pydantic import ValidationError

    cfg_path = _write_yaml(
        tmp_path / "unknown_layer.yaml",
        {
            "output_dir": str(tmp_path),
            "samples": [
                {
                    "name": "d0",
                    "layers": {
                        "nucleosomes": str(tmp_path / "n.csv"),
                        "4mC": str(tmp_path / "bogus.bed"),  # not a real layer
                    },
                },
            ],
        },
    )
    with pytest.raises(ValidationError, match="Unknown layer"):
        load_config(cfg_path)


def test_minimal_layers_just_nucleosomes(tmp_path):
    """A sample with ONLY nucleosomes should validate (5mC etc. are optional)."""
    cfg_path = _write_yaml(
        tmp_path / "minimal.yaml",
        {
            "output_dir": str(tmp_path),
            "samples": [
                {
                    "name": "d0",
                    "layers": {"nucleosomes": str(tmp_path / "n.csv")},
                },
            ],
        },
    )
    cfg = load_config(cfg_path)
    assert "nucleosomes" in cfg.samples[0].layers
    assert "5mC" not in cfg.samples[0].layers


def test_get_sample_by_name(tmp_path):
    """Convenience helper to look up sample by name."""
    cfg_path = _write_yaml(
        tmp_path / "two_samples.yaml",
        {
            "output_dir": str(tmp_path),
            "samples": [
                _minimal_sample_dict(tmp_path, "d0"),
                _minimal_sample_dict(tmp_path, "d4"),
            ],
        },
    )
    cfg = load_config(cfg_path)
    s = cfg.get_sample("d4")
    assert s.name == "d4"
    with pytest.raises(KeyError):
        cfg.get_sample("notarealsample")


def test_extra_field_in_sample_rejected(tmp_path):
    """Typos in sample keys should be caught (extra='forbid')."""
    from pydantic import ValidationError

    cfg_path = _write_yaml(
        tmp_path / "typo.yaml",
        {
            "output_dir": str(tmp_path),
            "samples": [
                {
                    "name": "d0",
                    "layers": {"nucleosomes": str(tmp_path / "n.csv")},
                    "typo_field": "oops",  # not in schema
                },
            ],
        },
    )
    with pytest.raises(ValidationError):
        load_config(cfg_path)
