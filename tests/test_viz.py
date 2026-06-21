"""Tests for public visualization helpers."""

from __future__ import annotations

import matplotlib
import pytest

matplotlib.use("Agg")

import matplotlib.pyplot as plt

from PACKAGE.db import FiberDatabase
from PACKAGE.viz import single_molecule_heatmap


def test_single_molecule_heatmap_renders_layers(tiny_db_path_with_layers, tmp_path):
    with FiberDatabase(tiny_db_path_with_layers) as db:
        fig = single_molecule_heatmap(db, "chr1", 100, 900, sample="d0")

    ax = fig.axes[0]
    assert ax.get_title() == "d0: chr1:100-900 (5mC)"
    assert len(ax.patches) == 5
    assert len(ax.collections) == 2  # one fiber baseline and one modification layer

    output = tmp_path / "single_molecule.png"
    fig.savefig(output)
    assert output.stat().st_size > 0
    plt.close(fig)


def test_single_molecule_heatmap_handles_empty_region(tiny_db_path_with_layers):
    with FiberDatabase(tiny_db_path_with_layers) as db:
        fig = single_molecule_heatmap(db, "chr1", 3000, 4000, sample="d0")

    assert fig.axes[0].texts[0].get_text() == "No overlapping fibers"
    plt.close(fig)


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"start": 10, "end": 10}, "0 <= start < end"),
        ({"start": 10, "end": 20, "layer": "bad"}, "layer must be one of"),
        ({"start": 10, "end": 20, "max_fibers": 0}, "at least 1"),
        ({"start": 10, "end": 20, "sample": "bad"}, "Unknown sample"),
    ],
)
def test_single_molecule_heatmap_validates_inputs(
    tiny_db_path_with_layers, kwargs, message
):
    with FiberDatabase(tiny_db_path_with_layers) as db:
        with pytest.raises(ValueError, match=message):
            single_molecule_heatmap(db, "chr1", **kwargs)
