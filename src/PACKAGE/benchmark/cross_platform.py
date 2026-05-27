"""Cross-platform concordance benchmark.

The central novel methodological claim of the paper: that the same biology, run on
ONT and PacBio with the same enzyme (Hia5), produces concordant results through
PACKAGE.

Computes, for matched serum-LIF samples:

  - Per-region methylation correlation (ONT vs PacBio)
  - Per-fiber nucleosome position agreement
  - MSP / accessibility-call agreement
  - Aggregate metaplots at known regulatory features

STATUS: Stub. This is the paper's headline cross-platform figure.
"""

from __future__ import annotations

from pathlib import Path


def run_cross_platform_benchmark(
    ont_db_path: Path,
    pacbio_db_path: Path,
    output_dir: Path,
    annotation: str = "CGI",
) -> None:
    """Compute concordance metrics and figures between matched ONT and PacBio databases.

    Outputs:
        - ``cross_platform_concordance.csv``: per-region correlation table
        - ``cross_platform_concordance.pdf``: figure (paper Fig 2b or main figure)
    """
    raise NotImplementedError
