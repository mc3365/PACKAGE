"""Visualization helpers.

Intentionally minimal. The single headline function is the per-molecule heatmap:
rows = individual fibers, columns = genomic position, color = modification state. This
is the visual signature of PACKAGE and likely produces a paper figure panel.

Submodules:
    single_molecule  Per-molecule heatmap and related helpers
"""

from PACKAGE.viz.annotation_matrix import (
    plot_annotation_heatmap,
    plot_annotation_metaplot,
    plot_centered_annotation_heatmap,
    plot_centered_metaplot,
)
from PACKAGE.viz.ecdf import plot_feature_ecdfs
from PACKAGE.viz.single_molecule import single_molecule_heatmap

__all__ = [
    "plot_annotation_heatmap",
    "plot_annotation_metaplot",
    "plot_centered_annotation_heatmap",
    "plot_centered_metaplot",
    "plot_feature_ecdfs",
    "single_molecule_heatmap",
]
