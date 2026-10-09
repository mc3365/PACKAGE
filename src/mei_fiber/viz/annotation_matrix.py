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
    interpolation: str = "hanning",
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

    # Match the original ComplexHeatmap color ramps closely.
    nuc_cmap = LinearSegmentedColormap.from_list(
        "package_nuc", ["white", "#63B8FF", "navy"]
    ).with_extremes(bad="#e5e5e5")
    met_cmap = LinearSegmentedColormap.from_list(
        "package_met", ["white", "lightcoral", "darkred"]
    ).with_extremes(bad="#e5e5e5")
    fig_height = max(6.0, min(16.0, nuc.shape[0] / 100))
    fig = plt.figure(figsize=(14.0, fig_height), constrained_layout=False)
    grid = fig.add_gridspec(
        1,
        5,
        width_ratios=[8.0, 1.2, 8.0, 1.2, 1.8],
        left=0.18,
        right=0.83,
        bottom=0.08,
        top=0.92,
        wspace=0.05,
    )
    axes = [
        fig.add_subplot(grid[0, 0]),
        fig.add_subplot(grid[0, 1]),
        fig.add_subplot(grid[0, 2]),
        fig.add_subplot(grid[0, 3]),
    ]
    color_axis = fig.add_subplot(grid[0, 4])
    color_axis.axis("off")
    im0 = axes[0].imshow(
        nuc,
        aspect="auto",
        interpolation=interpolation,
        cmap=nuc_cmap,
        vmin=0,
        vmax=1,
    )
    axes[1].barh(
        np.arange(len(nuc_mean)),
        nuc_mean,
        color="#4f94cd",
        edgecolor="none",
        height=1.0,
    )
    im1 = axes[2].imshow(
        met,
        aspect="auto",
        interpolation=interpolation,
        cmap=met_cmap,
        vmin=0,
        vmax=1,
    )
    axes[3].barh(
        np.arange(len(met_mean)),
        met_mean,
        color="#cd5555",
        edgecolor="none",
        height=1.0,
    )
    for ax in (axes[0], axes[2]):
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(True)
            spine.set_linewidth(0.8)
            spine.set_color("black")
    for ax in (axes[1], axes[3]):
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(True)
            spine.set_linewidth(0.8)
            spine.set_color("black")
    axes[1].set_ylim(axes[0].get_ylim())
    axes[3].set_ylim(axes[2].get_ylim())
    nuc_ticks = [0, nuc.shape[1] // 4, nuc.shape[1] // 2, 3 * nuc.shape[1] // 4, nuc.shape[1] - 1]
    met_ticks = [0, met.shape[1] // 4, met.shape[1] // 2, 3 * met.shape[1] // 4, met.shape[1] - 1]
    labels = [
        f"-{extension_bp / 1000:g} kb",
        f"-{extension_bp / 2000:g} kb",
        "TSS" if center_type == "TSS" else "Center",
        f"+{extension_bp / 2000:g} kb",
        f"+{extension_bp / 1000:g} kb",
    ]
    axes[0].set_xticks(nuc_ticks, labels=labels, rotation=90, fontsize=8)
    axes[2].set_xticks(met_ticks, labels=labels, rotation=90, fontsize=8)
    axes[0].axvline(nuc.shape[1] / 2, color="#4d4d4d", linestyle="--", linewidth=0.7)
    axes[2].axvline(met.shape[1] / 2, color="#4d4d4d", linestyle="--", linewidth=0.7)
    axes[0].set_title("Nucleosome Occupancy", fontsize=11, fontweight="bold")
    axes[2].set_title(met_label, fontsize=11, fontweight="bold")
    axes[1].set_xlabel("Nuc\nMean", fontsize=10)
    axes[3].set_xlabel("Met\nMean", fontsize=10)
    axes[1].set_xlim(0, 1)
    axes[3].set_xlim(0, 1)
    axes[1].set_xticks([0, 0.5, 1], labels=["0", "0.5", "1"], fontsize=8)
    axes[3].set_xticks([0, 0.5, 1], labels=["0", "0.5", "1"], fontsize=8)
    nuc_cax = color_axis.inset_axes([0.08, 0.54, 0.22, 0.10])
    met_cax = color_axis.inset_axes([0.08, 0.40, 0.22, 0.10])
    cbar0 = fig.colorbar(im0, cax=nuc_cax, ticks=[0, 0.5, 1])
    cbar1 = fig.colorbar(im1, cax=met_cax, ticks=[0, 0.5, 1])
    cbar0.ax.set_yticklabels(["0", "0.5", "1"])
    cbar1.ax.set_yticklabels(["0", "0.5", "1"])
    cbar0.ax.set_title("Nuc\nOccupancy", fontsize=9, loc="left")
    cbar1.ax.set_title("5mC", fontsize=9, loc="left")
    for cbar in (cbar0, cbar1):
        cbar.ax.tick_params(labelsize=7, length=2, pad=1)
    if title:
        fig.suptitle(title, fontsize=13, fontweight="bold", y=0.965)
    return fig


def plot_centered_metaplot(
    metaplot_csvs: dict[str, Path],
    *,
    annotation: str,
    extension_bp: int,
    center_label: str = "Center",
    sample_counts: dict[str, int] | None = None,
    signal_colors: dict[tuple[str, str], str] | None = None,
    title: str | None = None,
) -> Any:
    """Plot center-based metaplot profiles for one or more samples."""
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:  # pragma: no cover
        raise ImportError(
            'Metaplot figures require matplotlib. Install with `pip install ".[viz]"`.'
        ) from exc

    def _sample_role(sample: str) -> str | None:
        normalized = sample.lower().replace("-", "_")
        if normalized == "d0" or normalized.endswith("_d0"):
            return "d0"
        if normalized == "d4" or normalized.endswith("_d4"):
            return "d4"
        return None

    role_colors = {
        "d0": {"nuc": "#63B8FF", "met": "lightcoral"},
        "d4": {"nuc": "navy", "met": "darkred"},
    }
    n_samples = max(len(metaplot_csvs), 1)
    signal_colors = signal_colors or {}

    def _default_color(sample: str, sample_index: int, signal: str) -> str:
        role = _sample_role(sample)
        if role in role_colors:
            return role_colors[role][signal]
        if n_samples == 1:
            fraction = 0.55
        else:
            fraction = 0.35 + 0.55 * (sample_index / (n_samples - 1))
        cmap_name = "Blues" if signal == "nuc" else "Reds"
        return plt.get_cmap(cmap_name)(fraction)

    fig, ax = plt.subplots(figsize=(8, 5), constrained_layout=True)
    for sample_index, (sample, path) in enumerate(metaplot_csvs.items()):
        df = pd.read_csv(path)
        x = np.linspace(-extension_bp, extension_bp, len(df), endpoint=False)
        ax.plot(
            x,
            df["nuc_mean"],
            color=signal_colors.get((sample, "nuc"), _default_color(sample, sample_index, "nuc")),
            linewidth=1.0,
            alpha=0.8,
            label=f"{sample} Nucleosome",
        )
        ax.plot(
            x,
            df["met_mean"],
            color=signal_colors.get((sample, "met"), _default_color(sample, sample_index, "met")),
            linewidth=1.0,
            alpha=0.8,
            label=f"{sample} 5mC",
        )
    ext_kb = extension_bp / 1000
    ax.axvline(0, color="#4d4d4d", linestyle="--", linewidth=0.8)
    ax.annotate(
        center_label,
        xy=(0, 1),
        xycoords=("data", "axes fraction"),
        xytext=(4, -10),
        textcoords="offset points",
        ha="left",
        va="top",
        fontsize=9,
        fontstyle="italic",
        color="#4d4d4d",
    )
    ax.set_xticks(
        [-extension_bp, -extension_bp / 2, 0, extension_bp / 2, extension_bp],
        [f"-{ext_kb:g}kb", f"-{ext_kb / 2:g}kb", center_label, f"+{ext_kb / 2:g}kb", f"+{ext_kb:g}kb"],
    )
    ax.set_xlabel(f"Distance from {center_label} (bp)")
    ax.set_ylabel("Mean signal (fraction)")
    ax.set_ylim(bottom=0)
    ax.grid(True, color="#e5e5e5", linewidth=0.7)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.8)
        spine.set_color("black")
    ax.legend(
        title="Sample & Signal",
        frameon=False,
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        borderaxespad=0,
    )
    ax.set_title(
        title or f"{annotation} -- Center-Based Metaplot (V4)",
        fontsize=13,
        fontweight="bold",
        pad=24,
    )
    if sample_counts:
        subtitle = ", ".join(
            f"{sample} (n={sample_counts[sample]})"
            for sample in metaplot_csvs
            if sample in sample_counts
        )
        if subtitle:
            ax.text(
                0.5,
                1.005,
                subtitle,
                transform=ax.transAxes,
                ha="center",
                va="bottom",
                fontsize=10,
                color="#666666",
            )
    return fig
