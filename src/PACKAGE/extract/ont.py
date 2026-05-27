"""ONT extraction: orchestrates modkit + fibertools-rs ft extract.

The ONT pipeline:

1. **Dorado** has already produced a BAM with MM/ML modification tags
2. ``modkit extract`` reads those tags and outputs per-position 5mC/5hmC calls
3. ``ft extract`` (fibertools-rs) reads m6A and MSP calls, outputs BED12

This module wraps both with sensible defaults for ONT Hia5 Fiber-seq, validates inputs,
and writes a reproducibility manifest.

STATUS: Stub. To be filled in once user walks through their existing extraction workflow.
The full design will follow what they've been doing manually, just orchestrated.
"""

from __future__ import annotations

from PACKAGE.config import Config
from PACKAGE.utils import get_logger

log = get_logger(__name__)


def extract_samples(config: Config, samples: list[str] | None = None) -> None:
    """Run ONT extraction pipeline on one or more samples.

    For each sample, this will:
      1. Validate that the BAM has the expected modification tags (MM/ML)
      2. Run ``modkit extract`` with the configured parameters
      3. Run ``ft extract`` (fibertools-rs) with the configured parameters
      4. Apply any post-processing the user does (e.g. BED12 cleanup)
      5. Compute per-sample QC metrics
      6. Write the reproducibility manifest

    Args:
        config: Validated PACKAGE config.
        samples: Subset of sample names to extract; None means all in config.
    """
    samples = samples or config.sample_names
    log.info(f"ONT extraction: samples={samples}")
    # TODO(user-walkthrough): fill in based on the user's existing modkit + ft workflow
    raise NotImplementedError(
        "ONT extraction will be implemented after user walkthrough."
    )
