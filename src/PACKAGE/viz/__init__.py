"""Visualization helpers.

Intentionally minimal. The single headline function is the per-molecule heatmap:
rows = individual fibers, columns = genomic position, color = modification state. This
is the visual signature of PACKAGE and likely produces a paper figure panel.

Submodules:
    single_molecule  Per-molecule heatmap and related helpers
"""

from PACKAGE.viz.single_molecule import single_molecule_heatmap

__all__ = ["single_molecule_heatmap"]
