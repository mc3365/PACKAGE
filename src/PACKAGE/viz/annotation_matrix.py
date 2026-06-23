"""Plots for annotation-centered methylation/nucleosome matrices."""

from __future__ import annotations

import warnings
from pathlib import Path
from typing import Any

import numpy as np


def _load_matrix(path: Path) -> np.ndarray:
    if not Path(path).exists():
        raise FileNotFoundError(path)
    return np.load(path)["matrix"]


def plot_annotation_heatmap(
    nuc_body_npz: Path,
    met_body_npz: Path,
    *,
    title: str | None = None,
) -> Any:
    """Plot region-by-position nucleosome occupancy and 5mC body matrices."""
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:  # pragma: no cover
        raise ImportError(
            'Heatmap plots require matplotlib. Install with `pip install ".[viz]"`.'
        ) from exc

    nuc = _load_matrix(nuc_body_npz)
    met = _load_matrix(met_body_npz)
    fig, axes = plt.subplots(2, 1, figsize=(8.5, 5.6), sharex=True, constrained_layout=True)
    im0 = axes[0].imshow(nuc, aspect="auto", interpolation="nearest", cmap="Greys", vmin=0, vmax=1)
    im1 = axes[1].imshow(met, aspect="auto", interpolation="nearest", cmap="Reds", vmin=0, vmax=1)
    axes[0].set_ylabel("Regions")
    axes[1].set_ylabel("Regions")
    axes[1].set_xlabel("Normalized annotation body")
    axes[0].set_title("Nucleosome occupancy")
    axes[1].set_title("5mC methylation")
    fig.colorbar(im0, ax=axes[0], label="Occupancy")
    fig.colorbar(im1, ax=axes[1], label="Methylation")
    if title:
        fig.suptitle(title)
    return fig


def plot_annotation_metaplot(
    nuc_meta_npz: Path,
    met_meta_npz: Path,
    *,
    title: str | None = None,
) -> Any:
    """Plot average nucleosome occupancy and 5mC signal over flank-body-flank bins."""
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:  # pragma: no cover
        raise ImportError(
            'Metaplot figures require matplotlib. Install with `pip install ".[viz]"`.'
        ) from exc

    nuc = _load_matrix(nuc_meta_npz)
    met = _load_matrix(met_meta_npz)
    with warnings.catch_warnings(), np.errstate(invalid="ignore"):
        warnings.simplefilter("ignore", RuntimeWarning)
        nuc_mean = np.nanmean(nuc, axis=0)
        met_mean = np.nanmean(met, axis=0)
    x = np.arange(len(nuc_mean))

    fig, ax1 = plt.subplots(figsize=(7.5, 3.8), constrained_layout=True)
    ax1.plot(x, nuc_mean, color="#55585c", linewidth=1.8, label="Nucleosome")
    ax1.set_ylabel("Nucleosome occupancy")
    ax1.set_xlabel("Normalized flank-body-flank position")
    ax1.spines[["top", "right"]].set_visible(False)

    ax2 = ax1.twinx()
    ax2.plot(x, met_mean, color="#c43c39", linewidth=1.8, label="5mC")
    ax2.set_ylabel("5mC methylation")
    ax2.spines[["top"]].set_visible(False)

    for boundary in [100, 300]:
        ax1.axvline(boundary, color="#b8b8b8", linestyle="--", linewidth=0.9)
    if title:
        ax1.set_title(title)
    handles = ax1.get_lines() + ax2.get_lines()
    ax1.legend(handles=handles, labels=[h.get_label() for h in handles], frameon=False)
    return fig
