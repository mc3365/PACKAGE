"""Input validation and per-sample QC metrics.

This module is platform-agnostic — the same QC concepts apply to ONT and PacBio.

Two layers:

1. **Input validation**: cheap checks that should pass before extraction starts.
   Catches "did the user point us at the wrong BAM?" type errors early.

2. **QC metrics**: summary statistics computed during or after extraction. These
   become part of the reproducibility manifest and are useful for spotting bad
   samples (e.g. low per-fiber call rates suggesting a tag-encoding issue).

STATUS: Stub. To be filled in alongside the platform extractors.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any


def validate_bam(bam_path: Path, platform: str) -> None:
    """Validate that a BAM file looks like what we expect for the given platform.

    Checks:
      - File exists and is readable
      - Has expected modification tags (MM/ML for ONT; PacBio-specific for PacBio)
      - Is sorted and indexed
      - Reference name conventions match what the rest of the pipeline expects

    Raises:
        ValueError: with a clear message if validation fails.
    """
    raise NotImplementedError("Port from user's manual validation steps")


def compute_qc_metrics(sample_dir: Path) -> dict[str, Any]:
    """Compute per-sample QC metrics after extraction.

    Returns:
        Dict of metric name → value. Suggested metrics:
          - n_fibers: total number of fibers
          - mean_fiber_length: mean fiber length in bp
          - median_fiber_length
          - per_fiber_5mc_call_rate: fraction of CpGs called per fiber (median)
          - per_fiber_msp_count: median MSPs per fiber
          - chromosome_coverage: distribution of fibers across chromosomes
    """
    raise NotImplementedError("Define metrics during user walkthrough")
