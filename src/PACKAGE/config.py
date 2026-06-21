"""Typed configuration for PACKAGE.

Configuration is loaded from a YAML file rather than hardcoded in scripts. This gives:

  - One file that captures all parameters for a run (ships with paper supplements)
  - Type validation at load time, not at use time
  - Clear separation of "what to run" from "how to run it"

Example YAML
============

    output_dir: /path/to/output
    output_file: fiber_database.h5
    reference: /path/to/mm10/genome.fa

    samples:
      - name: esc_2i_d0
        bam: /path/to/d0.fiberseq.bam          # for the extract phase
        layers:
          nucleosomes: /path/to/d0_nuc_features.csv      # REQUIRED for build
          "5mC":       /path/to/d0_raw_mods.tsv.gz       # optional
          "5hmC":      /path/to/d0_raw_mods.tsv.gz       # SAME FILE as 5mC
          "6mA":       /path/to/d0_6ma.bed
          msp:         /path/to/d0_msp.bed

    annotations:
      master: /path/to/master_annotations_v3.bed

    parameters:
      methylation_threshold: 0.5
      min_msp_size: 10

    extraction:
      threads: 16
      overwrite: false
      keep_nucleosome_bed: true

Load with::

    from PACKAGE.config import load_config
    cfg = load_config("config.yaml")
"""

from __future__ import annotations

from pathlib import Path

import yaml
from pydantic import BaseModel, ConfigDict, Field, field_validator

# Layer names recognized by the config layer. Mirrors PACKAGE.db.schema.SUPPORTED_LAYERS
# but defined here to avoid a circular import (the config is used before the schema is
# loaded in some import orders).
_VALID_LAYER_NAMES: frozenset[str] = frozenset(
    {"nucleosomes", "5mC", "5hmC", "6mA", "msp"}
)


class SampleConfig(BaseModel):
    """One sample's input files for the build phase.

    Per-layer file paths are explicit so users with non-standard naming conventions
    are not constrained. Only ``nucleosomes`` is required (it defines which fibers
    exist on which chromosomes); other layers are optional and skipped if absent
    from ``layers``.
    """

    model_config = ConfigDict(extra="forbid")  # catch typos in YAML keys

    name: str = Field(description="Sample name (becomes the HDF5 top-level group name).")
    bam: Path | None = Field(
        default=None,
        description="Path to the fiberseq-annotated BAM (used by the extract phase).",
    )
    layers: dict[str, Path] = Field(
        description=(
            "Per-layer file paths. 'nucleosomes' is required; other layers "
            "(5mC, 5hmC, 6mA, msp) are optional and skipped if absent."
        ),
    )

    @field_validator("layers")
    @classmethod
    def _check_layer_names_and_nucleosomes(cls, v: dict[str, Path]) -> dict[str, Path]:
        unknown = set(v.keys()) - _VALID_LAYER_NAMES
        if unknown:
            raise ValueError(
                f"Unknown layer name(s) {sorted(unknown)} in sample.layers. "
                f"Valid: {sorted(_VALID_LAYER_NAMES)}."
            )
        if "nucleosomes" not in v:
            raise ValueError(
                "Sample must include a 'nucleosomes' layer (it defines which fibers "
                "exist on each chromosome). Even if you only care about other "
                "modifications, the nucleosomes file is required."
            )
        return v


class ParametersConfig(BaseModel):
    """Numerical parameters affecting database construction."""

    model_config = ConfigDict(extra="forbid")

    methylation_threshold: float = Field(
        default=0.5,
        ge=0.0,
        le=1.0,
        description="Per-CpG probability ≥ threshold → methylated (binary call).",
    )
    min_msp_size: int = Field(
        default=10,
        ge=0,
        description="Discard MSP blocks smaller than this (likely artifacts).",
    )


class ExtractionConfig(BaseModel):
    """External-tool settings for BAM-to-intermediate extraction."""

    model_config = ConfigDict(extra="forbid")

    threads: int = Field(default=16, ge=1)
    overwrite: bool = False
    keep_nucleosome_bed: bool = True
    validate_bam: bool = True
    write_manifest: bool = True
    modkit_executable: str = "modkit"
    ft_executable: str = "ft"
    samtools_executable: str = "samtools"


class Config(BaseModel):
    """Top-level PACKAGE configuration."""

    model_config = ConfigDict(extra="forbid")

    output_dir: Path = Field(description="Where to write the HDF5 database.")
    output_file: str = Field(
        default="fiber_database.h5",
        description="Filename for the output HDF5 database.",
    )
    reference: Path | None = Field(
        default=None,
        description="Path to the reference genome FASTA (used by extract phase).",
    )
    samples: list[SampleConfig] = Field(description="One entry per sample.")
    annotations: dict[str, Path] = Field(
        default_factory=dict,
        description=(
            "Mapping of annotation-set name → path to BED4 file. "
            "Currently 'master' is the only recognized name."
        ),
    )
    parameters: ParametersConfig = Field(default_factory=ParametersConfig)
    extraction: ExtractionConfig = Field(default_factory=ExtractionConfig)

    @property
    def output_path(self) -> Path:
        """Full path to the output HDF5 file."""
        return self.output_dir / self.output_file

    @property
    def sample_names(self) -> list[str]:
        """Convenience: list of sample names in declaration order."""
        return [s.name for s in self.samples]

    def get_sample(self, name: str) -> SampleConfig:
        """Look up a sample's config by name. Raises KeyError if not found."""
        for s in self.samples:
            if s.name == name:
                return s
        raise KeyError(f"Sample {name!r} not in config (have: {self.sample_names}).")


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
