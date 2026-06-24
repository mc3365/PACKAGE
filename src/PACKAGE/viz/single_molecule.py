"""Single-molecule visualization for genomic regions."""

from __future__ import annotations

from typing import Any

_SUPPORTED_LAYERS = {"5mC", "5hmC", "6mA"}
_SUPPORTED_TRACKS = {"full", "modification", "chromatin"}


def _get_fiber_span(
    db: Any, fiber_id: str, chrom: str, sample: str
) -> tuple[int, int] | None:
    """Return the genomic span for one fiber, or None if metadata is unavailable."""
    int_id = db._str_to_int(fiber_id, sample, chrom)
    if int_id is None:
        return None
    meta_path = f"{sample}/{chrom}/fiber_metadata"
    starts_path = f"{meta_path}/starts"
    if starts_path not in db.db:
        return None
    ids = db.db[f"{meta_path}/fiber_int_ids"][:]
    matches = ids == int(int_id)
    if not matches.any():
        return None
    idx = int(matches.nonzero()[0][0])
    start = int(db.db[starts_path][idx])
    end = int(db.db[f"{meta_path}/ends"][idx])
    return start, end


def single_molecule_heatmap(
    db: Any,
    chrom: str,
    start: int,
    end: int,
    sample: str | None = None,
    layer: str = "5mC",
    max_fibers: int = 100,
    show_unmodified: bool = True,
    tracks: str = "full",
    ax: Any = None,
) -> Any:
    """Render modification calls, nucleosomes, and MSPs for individual fibers.

    Each row represents one molecule. In ``tracks="full"`` mode, nucleosomes are
    gray blocks, MSPs are teal outlines, and modification calls are vertical marks.
    For 5mC and 5hmC, called modified bases are red and unmodified bases are light
    gray. 6mA calls are gold.
    The implementation draws sparse genomic features rather than allocating a dense
    base-by-fiber matrix.

    Args:
        db: An open FiberDatabase.
        chrom, start, end: Region to plot.
        sample: Sample name; if None, uses the first sample in the db.
        layer: Which modification to color by ('5mC', '6mA', '5hmC').
        max_fibers: Cap number of fibers shown to keep figures legible.
        show_unmodified: Show unmodified CpG calls for 5mC/5hmC layers.
        tracks: Track set to draw: ``full``, ``modification``, or ``chromatin``.
        ax: Optional existing matplotlib axes; if None, a new figure is created.

    Returns:
        The matplotlib ``Figure`` containing the plot.

    Raises:
        ValueError: If coordinates, sample, layer, or ``max_fibers`` are invalid.
        ImportError: If matplotlib is not installed.
    """
    try:
        import matplotlib.pyplot as plt
        from matplotlib.lines import Line2D
        from matplotlib.patches import Patch, Rectangle
    except ImportError as exc:  # pragma: no cover - depends on optional environment
        raise ImportError(
            'Visualization requires matplotlib. Install it with `pip install ".[viz]"`.'
        ) from exc

    if start < 0 or end <= start:
        raise ValueError("Region must satisfy 0 <= start < end")
    if layer not in _SUPPORTED_LAYERS:
        raise ValueError(f"layer must be one of {sorted(_SUPPORTED_LAYERS)}")
    if tracks not in _SUPPORTED_TRACKS:
        raise ValueError(f"tracks must be one of {sorted(_SUPPORTED_TRACKS)}")
    if max_fibers < 1:
        raise ValueError("max_fibers must be at least 1")

    sample = sample or db.samples[0]
    if sample not in db.samples:
        raise ValueError(f"Unknown sample {sample!r}; available samples: {db.samples}")

    fiber_ids = sorted(db.get_fibers_at(chrom, start, end, sample=sample))[:max_fibers]
    if ax is None:
        height = max(2.5, min(12.0, 1.5 + 0.22 * max(len(fiber_ids), 1)))
        fig, ax = plt.subplots(figsize=(12.5, height), constrained_layout=True)
    else:
        fig = ax.figure

    show_chromatin = tracks in {"full", "chromatin"}
    show_modification = tracks in {"full", "modification"}
    row_height = 0.62
    for row, fiber_id in enumerate(fiber_ids):
        fiber_span = _get_fiber_span(db, fiber_id, chrom, sample)
        if fiber_span is None:
            fiber_left, fiber_right = start, end
        else:
            fiber_left = max(fiber_span[0], start)
            fiber_right = min(fiber_span[1], end)
        if fiber_right > fiber_left:
            ax.hlines(
                row,
                fiber_left,
                fiber_right,
                color="#111111",
                linewidth=0.45,
                zorder=0,
            )

        if show_chromatin:
            nucleosomes = db.get_nucleosomes(fiber_id, chrom, sample=sample)
            for feature_start, feature_end in zip(
                nucleosomes.get("starts", []), nucleosomes.get("ends", []), strict=False
            ):
                left = max(int(feature_start), start)
                right = min(int(feature_end), end)
                if right > left:
                    ax.add_patch(
                        Rectangle(
                            (left, row - row_height / 2),
                            right - left,
                            row_height,
                            facecolor="#55585c",
                            edgecolor="none",
                            alpha=0.9,
                            zorder=1,
                        )
                    )

            msps = db.get_msp(fiber_id, chrom, sample=sample)
            for feature_start, feature_end in zip(
                msps.get("starts", []), msps.get("ends", []), strict=False
            ):
                left = max(int(feature_start), start)
                right = min(int(feature_end), end)
                if right > left:
                    ax.add_patch(
                        Rectangle(
                            (left, row - row_height / 2),
                            right - left,
                            row_height,
                            facecolor="none",
                            edgecolor="#168a8a",
                            linewidth=0.85,
                            alpha=0.9,
                            zorder=2,
                        )
                    )

        if not show_modification:
            continue
        modifications = db.get_methylation(
            fiber_id, chrom, mod_type=layer, sample=sample
        )
        positions = modifications.get("positions", [])
        in_region = [i for i, position in enumerate(positions) if start <= position < end]
        if not in_region:
            continue

        plot_positions = [positions[i] for i in in_region]
        if layer == "6mA" or "is_methylated" not in modifications:
            ax.scatter(
                plot_positions,
                [row] * len(in_region),
                marker="|",
                s=46,
                linewidths=1.2,
                c="#e3a018",
                alpha=0.95,
                zorder=3,
            )
        else:
            calls = modifications["is_methylated"]
            modified_positions = [
                positions[i] for i in in_region if bool(calls[i])
            ]
            unmodified_positions = [
                positions[i] for i in in_region if not bool(calls[i])
            ]
            if show_unmodified and unmodified_positions:
                ax.scatter(
                    unmodified_positions,
                    [row] * len(unmodified_positions),
                    marker="|",
                    s=30,
                    linewidths=0.8,
                    c="#b7bdc3",
                    alpha=0.4,
                    zorder=3,
                )
            if modified_positions:
                ax.scatter(
                    modified_positions,
                    [row] * len(modified_positions),
                    marker="|",
                    s=48,
                    linewidths=1.2,
                    c="#c43c39",
                    alpha=0.95,
                    zorder=4,
                )

    ax.set_xlim(start, end)
    ax.set_ylim(-0.75, max(len(fiber_ids) - 0.25, 0.75))
    ax.invert_yaxis()
    ax.set_xlabel(f"{chrom} position (bp)")
    ax.set_ylabel("Individual fibers")
    if tracks == "chromatin":
        title_label = "chromatin"
    elif tracks == "modification":
        title_label = f"{layer}, modification"
    else:
        title_label = layer
    ax.set_title(f"{sample}: {chrom}:{start:,}-{end:,} ({title_label})")
    ax.ticklabel_format(axis="x", style="plain", useOffset=False)
    ax.set_yticks([])
    ax.spines[["top", "right", "left"]].set_visible(False)

    if not fiber_ids:
        ax.text(
            0.5,
            0.5,
            "No overlapping fibers",
            transform=ax.transAxes,
            ha="center",
            va="center",
            color="#55585c",
        )

    handles = []
    if show_chromatin:
        handles.extend(
            [
                Patch(facecolor="#55585c", label="Nucleosome"),
                Patch(facecolor="none", edgecolor="#168a8a", label="MSP"),
            ]
        )
    if show_modification:
        modification_handles = (
            [Line2D([], [], color="#e3a018", marker="|", linestyle="None", label="6mA")]
            if layer == "6mA"
            else [Line2D([], [], color="#c43c39", marker="|", linestyle="None", label=layer)]
        )
        if layer != "6mA" and show_unmodified:
            modification_handles.append(
                Line2D([], [], color="#b7bdc3", marker="|", linestyle="None", label="Unmodified")
            )
        handles.extend(modification_handles)
    if handles:
        ax.legend(
            handles=handles,
            loc="upper left",
            bbox_to_anchor=(1.005, 1.0),
            borderaxespad=0,
            frameon=False,
        )
    return fig
