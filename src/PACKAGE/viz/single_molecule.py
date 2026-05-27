"""Per-molecule visualization: heatmap of individual fibers in a region.

The flagship visualization for PACKAGE. Shows what bulk methods cannot: each row is
one molecule, columns are positions in a genomic region, and color encodes the
modification state at that position on that molecule. Nucleosome positions and MSPs
can be overlaid as horizontal bars.

STATUS: Stub. Implementation deferred until the database query layer is solid.
"""

from __future__ import annotations

from typing import Any


def single_molecule_heatmap(
    db: Any,                            # PACKAGE.db.FiberDatabase
    chrom: str,
    start: int,
    end: int,
    sample: str | None = None,
    layer: str = "5mC",
    max_fibers: int = 100,
    ax: Any = None,                     # matplotlib.axes.Axes | None
) -> Any:                               # matplotlib.figure.Figure
    """Render a per-molecule heatmap for one genomic region.

    Args:
        db: An open FiberDatabase.
        chrom, start, end: Region to plot.
        sample: Sample name; if None, uses the first sample in the db.
        layer: Which modification to color by ('5mC', '6mA', '5hmC').
        max_fibers: Cap number of fibers shown to keep figures legible.
        ax: Optional existing matplotlib axes; if None, a new figure is created.

    Returns:
        The matplotlib Figure.
    """
    raise NotImplementedError("Single-molecule heatmap to be implemented")
