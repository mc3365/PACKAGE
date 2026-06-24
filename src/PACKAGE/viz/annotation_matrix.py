"""Plots for annotation-centered methylation/nucleosome matrices."""

from __future__ import annotations

import warnings
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


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


def plot_centered_annotation_heatmap(
    nuc_npz: Path,
    met_npz: Path,
    *,
    met_label: str = "5mC",
    title: str | None = None,
) -> Any:
    """Plot paired center-based nucleosome and 5mC heatmaps."""
    try:
        import matplotlib.pyplot as plt
        from matplotlib.colors import LinearSegmentedColormap
    except ImportError as exc:  # pragma: no cover
        raise ImportError(
            'Heatmap plots require matplotlib. Install with `pip install ".[viz]"`.'
        ) from exc

    nuc_data = np.load(nuc_npz, allow_pickle=True)
    met_data = np.load(met_npz, allow_pickle=True)
    nuc = nuc_data["matrix"]
    met = met_data["matrix"]
    extension_bp = int(nuc_data["extension_bp"])
    center_type = str(nuc_data["center_type"])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        order = np.argsort(-np.nanmean(nuc, axis=1))
        nuc_mean = np.nanmean(nuc[order], axis=1)
        met_mean = np.nanmean(met[order], axis=1)
    nuc = nuc[order]
    met = met[order]

    nuc_cmap = LinearSegmentedColormap.from_list("package_nuc", ["white", "lightskyblue", "navy"])
    met_cmap = LinearSegmentedColormap.from_list("package_met", ["white", "lightcoral", "darkred"])
    fig, axes = plt.subplots(
        1,
        4,
        figsize=(11.5, max(4.5, min(12.0, nuc.shape[0] / 80))),
        width_ratios=[5.0, 0.8, 5.0, 0.8],
        constrained_layout=True,
    )
    im0 = axes[0].imshow(nuc, aspect="auto", interpolation="nearest", cmap=nuc_cmap, vmin=0, vmax=1)
    axes[1].barh(np.arange(len(nuc_mean)), nuc_mean, color="steelblue")
    im1 = axes[2].imshow(met, aspect="auto", interpolation="nearest", cmap=met_cmap, vmin=0, vmax=1)
    axes[3].barh(np.arange(len(met_mean)), met_mean, color="indianred")
    for ax in axes:
        ax.invert_yaxis()
        ax.set_yticks([])
        ax.spines[["top", "right", "left"]].set_visible(False)
    nuc_ticks = [0, nuc.shape[1] // 4, nuc.shape[1] // 2, 3 * nuc.shape[1] // 4, nuc.shape[1] - 1]
    met_ticks = [0, met.shape[1] // 4, met.shape[1] // 2, 3 * met.shape[1] // 4, met.shape[1] - 1]
    labels = [
        f"-{extension_bp / 1000:g} kb",
        f"-{extension_bp / 2000:g} kb",
        "TSS" if center_type == "TSS" else "Center",
        f"+{extension_bp / 2000:g} kb",
        f"+{extension_bp / 1000:g} kb",
    ]
    axes[0].set_xticks(nuc_ticks, labels=labels)
    axes[2].set_xticks(met_ticks, labels=labels)
    axes[0].axvline(nuc.shape[1] / 2, color="#4d4d4d", linestyle="--", linewidth=0.7)
    axes[2].axvline(met.shape[1] / 2, color="#4d4d4d", linestyle="--", linewidth=0.7)
    axes[0].set_title("Nucleosome Occupancy")
    axes[2].set_title(met_label)
    axes[1].set_title("Nuc\nMean")
    axes[3].set_title("Met\nMean")
    axes[1].set_xlim(0, 1)
    axes[3].set_xlim(0, 1)
    fig.colorbar(im0, ax=axes[0], label="Occupancy")
    fig.colorbar(im1, ax=axes[2], label="5mC")
    if title:
        fig.suptitle(title)
    return fig


def plot_centered_metaplot(
    metaplot_csvs: dict[str, Path],
    *,
    annotation: str,
    extension_bp: int,
    center_label: str = "Center",
    title: str | None = None,
) -> Any:
    """Plot center-based metaplot profiles with four-color d0/d4 styling."""
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:  # pragma: no cover
        raise ImportError(
            'Metaplot figures require matplotlib. Install with `pip install ".[viz]"`.'
        ) from exc

    colors = {
        "d0_nuc": "lightskyblue",
        "d4_nuc": "navy",
        "d0_met": "lightcoral",
        "d4_met": "darkred",
    }
    fig, ax = plt.subplots(figsize=(8, 5), constrained_layout=True)
    for sample, path in metaplot_csvs.items():
        df = pd.read_csv(path)
        x = np.linspace(-extension_bp, extension_bp, len(df), endpoint=False)
        ax.plot(
            x,
            df["nuc_mean"],
            color=colors.get(f"{sample}_nuc", "steelblue"),
            linewidth=1.4,
            label=f"{sample} Nucleosome",
        )
        ax.plot(
            x,
            df["met_mean"],
            color=colors.get(f"{sample}_met", "darkred"),
            linewidth=1.4,
            label=f"{sample} 5mC",
        )
    ext_kb = extension_bp / 1000
    ax.axvline(0, color="#4d4d4d", linestyle="--", linewidth=0.8)
    ax.set_xticks(
        [-extension_bp, -extension_bp / 2, 0, extension_bp / 2, extension_bp],
        [f"-{ext_kb:g}kb", f"-{ext_kb / 2:g}kb", center_label, f"+{ext_kb / 2:g}kb", f"+{ext_kb:g}kb"],
    )
    ax.set_xlabel(f"Distance from {center_label} (bp)")
    ax.set_ylabel("Mean signal (fraction)")
    ax.set_ylim(bottom=0)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(title="Sample & Signal", frameon=False, loc="best")
    ax.set_title(title or f"{annotation} center-based metaplot")
    return fig
