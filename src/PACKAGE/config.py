"""Typed configuration for PACKAGE.

Configuration is loaded from a YAML file rather than hardcoded in scripts. This gives:

  - One file that captures all parameters for a run (ships with paper supplements)
  - Type validation at load time, not at use time
  - Clear separation of "what to run" from "how to run it"

Example YAML:

    # config.yaml
    output_dir: /vast/palmer/scratch/zsmith/mc3365/long_read/h5
    output_file: fiber_database.h5
    samples:
      - name: d0
        input_dir: /path/to/d0
        ft_extract: d0.fibers.bed.gz
        modkit_extract: d0.mods.tsv.gz
      - name: d4
        input_dir: /path/to/d4
        ft_extract: d4.fibers.bed.gz
        modkit_extract: d4.mods.tsv.gz
    annotations:
      master_v3: /path/to/master_annotations_v3.bed
    parameters:
      methylation_threshold: 0.5
      min_msp_size: 10

Load with:

    from PACKAGE.config import load_config
    cfg = load_config("config.yaml")
"""

from __future__ import annotations

from pathlib import Path

import yaml
from pydantic import BaseModel, Field


class SampleConfig(BaseModel):
    """One sample's input file locations."""

    name: str = Field(description="Sample name used as HDF5 group key (e.g. 'd0', 'd4').")
    input_dir: Path = Field(description="Directory containing this sample's input files.")
    ft_extract: str = Field(description="Filename of ft-extract output (BED12).")
    modkit_extract: str | None = Field(
        default=None,
        description="Filename of modkit extract output (ONT only).",
    )
    # Add additional per-sample fields here as needed (e.g. pb_cpg_tools output for PacBio).


class ParametersConfig(BaseModel):
    """Numerical parameters affecting database construction."""

    methylation_threshold: float = Field(
        default=0.5,
        ge=0.0,
        le=1.0,
        description="Probability ≥ threshold → methylated (binary call).",
    )
    min_msp_size: int = Field(
        default=10,
        ge=0,
        description="Discard MSP blocks smaller than this (likely artifacts).",
    )


class Config(BaseModel):
    """Top-level PACKAGE configuration."""

    output_dir: Path = Field(description="Where to write the HDF5 database.")
    output_file: str = Field(
        default="fiber_database.h5",
        description="Filename for the output HDF5 database.",
    )
    samples: list[SampleConfig] = Field(description="One entry per sample.")
    annotations: dict[str, Path] = Field(
        default_factory=dict,
        description="Mapping of annotation-set name → path to BED/TSV file.",
    )
    parameters: ParametersConfig = Field(default_factory=ParametersConfig)

    @property
    def output_path(self) -> Path:
        """Full path to the output HDF5 file."""
        return self.output_dir / self.output_file

    @property
    def sample_names(self) -> list[str]:
        """Convenience: list of sample names in order."""
        return [s.name for s in self.samples]


def load_config(path: str | Path) -> Config:
    """Load and validate a PACKAGE YAML config file.

    Raises:
        FileNotFoundError: if the config file does not exist.
        pydantic.ValidationError: if the config is malformed.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Config file not found: {path}")
    with path.open() as f:
        raw = yaml.safe_load(f)
    return Config.model_validate(raw)
