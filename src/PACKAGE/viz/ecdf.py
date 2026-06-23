"""ECDF plots for per-fiber feature fractions."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

_COVERAGE_COLUMNS = {
    "frac_5mC_coverage": "5mC",
    "frac_5hmC_coverage": "5hmC",
    "frac_6mA_coverage": "6mA",
    "frac_msp_coverage": "MSP",
    "frac_nuc_coverage": "Nucleosome",
}
_BASE_COLUMNS = {
    "frac_5mC_among_cpg": "5mC among CpG",
    "frac_5hmC_among_cpg": "5hmC among CpG",
    "frac_6mA_among_A": "6mA density proxy",
}


def _ecdf_xy(values: pd.Series, x_min: float = 0.0, x_max: float = 1.0) -> tuple[np.ndarray, np.ndarray]:
    """Return x/y coordinates for an ECDF line that spans the full plot range."""
    x = values.dropna().sort_values().to_numpy(dtype=float)
    if len(x) == 0:
        return np.array([], dtype=float), np.array([], dtype=float)
    y = np.arange(1, len(x) + 1, dtype=float) / len(x)
    x = np.concatenate([[x_min], x, [x_max]])
    y = np.concatenate([[0.0], y, [1.0]])
    return x, y


def _plot_ecdf_set(
    df: pd.DataFrame,
    columns: dict[str, str],
    *,
    title: str,
    xlabel: str,
    ncols: int,
) -> Any:
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:  # pragma: no cover
        raise ImportError(
            'ECDF plots require matplotlib. Install with `pip install ".[viz]"`.'
        ) from exc

    samples = list(dict.fromkeys(df["sample"].astype(str)))
    colors = ["#3b6fb6", "#c74f46", "#2f8f6b", "#8d63b8", "#8a7a2f"]
    nrows = (len(columns) + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(3.4 * ncols, 3.0 * nrows), squeeze=False)
    for ax, (column, label) in zip(axes.ravel(), columns.items(), strict=False):
        for sample, color in zip(samples, colors, strict=False):
            values = df.loc[df["sample"].astype(str) == sample, column].dropna().sort_values()
            if values.empty:
                continue
            x, y = _ecdf_xy(values)
            ax.step(x, y, where="post", color=color, linewidth=1.2, label=sample)
        ax.set_title(label)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.spines[["top", "right"]].set_visible(False)
    for ax in axes.ravel()[len(columns):]:
        ax.axis("off")
    axes[0, 0].set_ylabel("ECDF")
    for ax in axes[-1, :]:
        ax.set_xlabel(xlabel)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.suptitle(title, y=0.99, fontsize=13)
    if handles:
        fig.legend(
            handles,
            labels,
            loc="upper center",
            bbox_to_anchor=(0.5, 0.94),
            ncols=min(len(samples), 4),
            frameon=False,
        )
    fig.tight_layout(rect=(0, 0, 1, 0.86))
    return fig


def plot_feature_ecdfs(feature_csv: Path) -> tuple[Any, Any]:
    """Plot coverage-fraction and base-specific ECDF panels from a feature CSV."""
    df = pd.read_csv(feature_csv)
    missing = [
        col for col in ["sample", *_COVERAGE_COLUMNS, *_BASE_COLUMNS]
        if col not in df.columns
    ]
    if missing:
        raise ValueError(f"Feature CSV is missing columns: {missing}")
    coverage_fig = _plot_ecdf_set(
        df,
        _COVERAGE_COLUMNS,
        title="Per-fiber feature coverage",
        xlabel="Fraction of fiber length",
        ncols=3,
    )
    base_fig = _plot_ecdf_set(
        df,
        _BASE_COLUMNS,
        title="Per-fiber base-call fractions",
        xlabel="Fraction",
        ncols=3,
    )
    return coverage_fig, base_fig
