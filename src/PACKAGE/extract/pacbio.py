"""PacBio extraction: orchestrates pb-CpG-tools or fibertools (5mC) + fibertools-rs ft extract.

The PacBio pipeline differs from ONT in three places:

1. **5mC calls** come from pb-CpG-tools or fibertools-rs native 5mC mode (not modkit)
2. **Tag conventions** in the BAM are PacBio-specific (kinetics-based mod calling)
3. **m6A and MSP** are computed by fibertools-rs trained models specific to PacBio data

The schema-level data layout is identical to ONT — what differs is which mod types are
populated and the upstream tool used.

STATUS: Stub. To be filled in once user walks through their existing extraction workflow,
including how the lab's existing serum-LIF PacBio data was processed.
"""

from __future__ import annotations

from PACKAGE.config import Config
from PACKAGE.utils import get_logger

log = get_logger(__name__)


def extract_samples(config: Config, samples: list[str] | None = None) -> None:
    """Run PacBio extraction pipeline on one or more samples.

    Args:
        config: Validated PACKAGE config.
        samples: Subset of sample names to extract; None means all in config.
    """
    samples = samples or config.sample_names
    log.info(f"PacBio extraction: samples={samples}")
    # TODO(user-walkthrough): fill in based on the lab's existing PacBio workflow
    raise NotImplementedError(
        "PacBio extraction will be implemented after user walkthrough."
    )
