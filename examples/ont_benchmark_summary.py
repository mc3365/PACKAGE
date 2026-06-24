#!/usr/bin/env python3
"""Create summary figures from PACKAGE ONT benchmark outputs."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any


GB = 1024**3
LAYER_ORDER = ["fibers", "5mC", "5hmC", "6mA", "msp", "nucleosomes"]
LAYER_COLORS = {
    "fibers": "#6b7280",
    "5mC": "#c43c39",
    "5hmC": "#f08a84",
    "6mA": "#e3a018",
    "msp": "#168a8a",
    "nucleosomes": "#55585c",
}


def _require_matplotlib() -> Any:
    try:
        import matplotlib.pyplot as plt
        from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
    except ImportError as exc:  # pragma: no cover
        raise ImportError(
            'Benchmark summary plots require matplotlib. Install with `pip install ".[viz]"`.'
        ) from exc
    return plt, FancyArrowPatch, FancyBboxPatch


def _gb(size_bytes: int | float | None) -> float:
    return 0.0 if size_bytes is None else float(size_bytes) / GB


def plot_storage_summary(summary: dict[str, Any], outdir: Path, bam_paths: list[Path]) -> None:
    plt, _, _ = _require_matplotlib()

    labels = []
    sizes = []
    if bam_paths:
        labels.append("Input BAM")
        sizes.append(sum(path.stat().st_size for path in bam_paths if path.exists()))
    labels.extend(["Extracted\nintermediates", "PACKAGE\nHDF5", "Spatial\nindex"])
    sizes.extend(
        [
            int(summary.get("intermediate_size_bytes", 0)),
            int(summary.get("database_size_bytes", 0)),
            int(summary.get("spatial_index_size_bytes", 0)),
        ]
    )

    fig, ax = plt.subplots(figsize=(6.7, 4.2), constrained_layout=True)
    colors = ["#9ca3af", "#4c78a8", "#59a14f", "#b07aa1"][-len(labels):]
    bars = ax.bar(labels, [_gb(size) for size in sizes], color=colors)
    ax.set_ylabel("Size (GB)")
    ax.set_title("Storage footprint")
    ax.spines[["top", "right"]].set_visible(False)
    for bar, size in zip(bars, sizes, strict=False):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height(),
            f"{_gb(size):.1f} GB",
            ha="center",
            va="bottom",
            fontsize=9,
        )

    intermediate = int(summary.get("intermediate_size_bytes", 0))
    database = int(summary.get("database_size_bytes", 0))
    if intermediate and database:
        ratio = intermediate / database
        reduction = 100 * (1 - database / intermediate)
        ax.text(
            0.98,
            0.93,
            f"{ratio:.1f}x smaller\n{reduction:.1f}% reduction",
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=10,
        )

    fig.savefig(outdir / "benchmark_storage_summary.png", dpi=220)
    fig.savefig(outdir / "benchmark_storage_summary.pdf")
    plt.close(fig)


def plot_record_counts(summary: dict[str, Any], outdir: Path) -> None:
    plt, _, _ = _require_matplotlib()

    samples = sorted(summary.get("samples", {}))
    fig, ax = plt.subplots(figsize=(7.2, 4.4), constrained_layout=True)
    bottoms = [0.0] * len(samples)
    for layer in LAYER_ORDER:
        values = [
            summary["samples"][sample].get(layer, 0) / 1e9
            for sample in samples
        ]
        if not any(values):
            continue
        ax.bar(samples, values, bottom=bottoms, label=layer, color=LAYER_COLORS.get(layer))
        bottoms = [old + value for old, value in zip(bottoms, values, strict=False)]

    ax.set_ylabel("Records (billions)")
    ax.set_title("Database scale by sample and layer")
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, ncols=3, loc="upper center", bbox_to_anchor=(0.5, -0.12))
    fig.savefig(outdir / "benchmark_record_counts.png", dpi=220)
    fig.savefig(outdir / "benchmark_record_counts.pdf")
    plt.close(fig)


def plot_random_query_results(results_csv: Path, outdir: Path, outlier_fibers: int) -> None:
    plt, _, _ = _require_matplotlib()
    import pandas as pd

    df = pd.read_csv(results_csv)
    if df.empty:
        return
    df["window_kb"] = df["region_size_bp"] / 1000
    df["elapsed_ms"] = df["elapsed_seconds"] * 1000
    df["is_outlier"] = df["n_fibers"] >= outlier_fibers

    fig, ax = plt.subplots(figsize=(7.2, 4.8), constrained_layout=True)
    typical = df[~df["is_outlier"]]
    outliers = df[df["is_outlier"]]
    sc = ax.scatter(
        typical["window_kb"],
        typical["elapsed_ms"],
        c=typical["n_fibers"],
        s=45,
        cmap="viridis",
        alpha=0.85,
        edgecolor="none",
        label="typical regions",
    )
    if not outliers.empty:
        ax.scatter(
            outliers["window_kb"],
            outliers["elapsed_ms"],
            c=outliers["n_fibers"],
            s=90,
            cmap="viridis",
            marker="D",
            edgecolor="#b00020",
            linewidth=1.3,
            label=f">= {outlier_fibers:,} fibers",
        )
        for _, row in outliers.iterrows():
            ax.annotate(
                f"{row['chrom']}\n{int(row['n_fibers']):,} fibers",
                (row["window_kb"], row["elapsed_ms"]),
                xytext=(8, 6),
                textcoords="offset points",
                fontsize=8,
            )

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Query window (kb)")
    ax.set_ylabel("Elapsed time (ms)")
    ax.set_title("Random-window coordinate queries")
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, loc="upper left")
    cbar = fig.colorbar(sc, ax=ax)
    cbar.set_label("Returned fibers")
    fig.savefig(outdir / "benchmark_random_query_scatter.png", dpi=220)
    fig.savefig(outdir / "benchmark_random_query_scatter.pdf")
    plt.close(fig)

    trimmed = df[~df["is_outlier"]]
    if not trimmed.empty:
        summary = (
            trimmed.groupby(["sample", "mode", "region_size_bp"], as_index=False)
            .agg(
                n_queries=("elapsed_seconds", "size"),
                median_seconds=("elapsed_seconds", "median"),
                q25_seconds=("elapsed_seconds", lambda values: values.quantile(0.25)),
                q75_seconds=("elapsed_seconds", lambda values: values.quantile(0.75)),
                max_n_fibers=("n_fibers", "max"),
            )
        )
        summary.to_csv(outdir / "query_speed_summary_without_outliers.csv", index=False)


def plot_annotation_benchmark(annotation_summary_csv: Path, outdir: Path) -> None:
    plt, _, _ = _require_matplotlib()
    import pandas as pd

    df = pd.read_csv(annotation_summary_csv)
    if df.empty:
        return
    df = df.sort_values("median_seconds")

    med = df["median_seconds"].astype(float)
    q25 = df["q25_seconds"].astype(float)
    q75 = df["q75_seconds"].astype(float)
    lower = med - q25
    upper = q75 - med

    fig, ax = plt.subplots(figsize=(7.0, 4.2), constrained_layout=True)
    bars = ax.bar(df["annotation"], med, color="#4c78a8")
    ax.errorbar(
        df["annotation"],
        med,
        yerr=[lower, upper],
        fmt="none",
        color="#1f2937",
        capsize=3,
        linewidth=1.0,
    )
    ax.set_ylabel("Median query time (s)")
    ax.set_title("Annotation query benchmark")
    ax.spines[["top", "right"]].set_visible(False)
    for bar, rows in zip(bars, df["n_rows_max"], strict=False):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height(),
            f"{int(rows):,} rows",
            ha="center",
            va="bottom",
            fontsize=9,
        )
    fig.savefig(outdir / "benchmark_annotation_queries.png", dpi=220)
    fig.savefig(outdir / "benchmark_annotation_queries.pdf")
    plt.close(fig)


def plot_workflow(outdir: Path) -> None:
    plt, FancyArrowPatch, FancyBboxPatch = _require_matplotlib()

    fig, ax = plt.subplots(figsize=(11, 4.2), constrained_layout=True)
    ax.set_axis_off()
    boxes = [
        ("ONT fiberseq\nBAM", 0.04, 0.55, "#dbeafe"),
        ("PACKAGE extract\nmodkit + ft", 0.24, 0.55, "#e0f2fe"),
        ("Intermediate files\n5mC/5hmC, 6mA,\nMSP, nucleosomes", 0.44, 0.55, "#fef3c7"),
        ("PACKAGE build\nHDF5 + indices", 0.64, 0.55, "#dcfce7"),
        ("Query + figures\nregions, annotations,\nsingle molecules", 0.84, 0.55, "#f3e8ff"),
    ]
    for text, x, y, color in boxes:
        patch = FancyBboxPatch(
            (x, y),
            0.13,
            0.26,
            boxstyle="round,pad=0.02,rounding_size=0.02",
            facecolor=color,
            edgecolor="#374151",
            linewidth=1.0,
            transform=ax.transAxes,
        )
        ax.add_patch(patch)
        ax.text(
            x + 0.065,
            y + 0.13,
            text,
            ha="center",
            va="center",
            fontsize=10,
            transform=ax.transAxes,
        )

    for x0, x1 in [(0.17, 0.24), (0.37, 0.44), (0.57, 0.64), (0.77, 0.84)]:
        arrow = FancyArrowPatch(
            (x0, 0.68),
            (x1, 0.68),
            arrowstyle="-|>",
            mutation_scale=14,
            linewidth=1.2,
            color="#374151",
            transform=ax.transAxes,
        )
        ax.add_patch(arrow)

    ax.text(
        0.5,
        0.22,
        (
            "Current validated stage: ONT extraction, HDF5 packaging, "
            "spatial-indexed queries, benchmarks, and visualization."
        ),
        ha="center",
        va="center",
        fontsize=11,
        transform=ax.transAxes,
    )
    fig.savefig(outdir / "package_ont_workflow.png", dpi=220)
    fig.savefig(outdir / "package_ont_workflow.pdf")
    plt.close(fig)


def write_text_summary(summary: dict[str, Any], outdir: Path, bam_paths: list[Path]) -> None:
    database = int(summary.get("database_size_bytes", 0))
    intermediates = int(summary.get("intermediate_size_bytes", 0))
    index = int(summary.get("spatial_index_size_bytes", 0))
    bam_size = sum(path.stat().st_size for path in bam_paths if path.exists())
    lines = [
        "# PACKAGE ONT Benchmark Summary",
        "",
        (
            f"- Input BAM size: {_gb(bam_size):.2f} GB"
            if bam_paths
            else "- Input BAM size: not provided"
        ),
        f"- Extracted intermediate size: {_gb(intermediates):.2f} GB",
        f"- HDF5 database size: {_gb(database):.2f} GB",
        f"- Spatial index size: {_gb(index):.2f} GB",
    ]
    if intermediates and database:
        lines.append(f"- Intermediate/HDF5 ratio: {intermediates / database:.2f}x")
        lines.append(
            f"- Storage reduction from intermediates: "
            f"{100 * (1 - database / intermediates):.2f}%"
        )
    if database and index:
        lines.append(f"- Spatial-index overhead relative to HDF5: {100 * index / database:.2f}%")
    total_records = 0
    for sample_totals in summary.get("samples", {}).values():
        total_records += sum(int(value) for value in sample_totals.values())
    if summary.get("annotation_records"):
        total_records += int(summary["annotation_records"])
    lines.append(f"- Stored fiber/layer/annotation records: {total_records:,}")
    (outdir / "benchmark_summary.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--storage-summary", type=Path, required=True)
    parser.add_argument("--query-results", type=Path)
    parser.add_argument("--annotation-summary", type=Path)
    parser.add_argument("--bam", type=Path, action="append", default=[])
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--outlier-fibers", type=int, default=10_000)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    summary = json.loads(args.storage_summary.read_text())
    plot_storage_summary(summary, args.outdir, args.bam)
    plot_record_counts(summary, args.outdir)
    if args.query_results:
        plot_random_query_results(args.query_results, args.outdir, args.outlier_fibers)
    if args.annotation_summary:
        plot_annotation_benchmark(args.annotation_summary, args.outdir)
    plot_workflow(args.outdir)
    write_text_summary(summary, args.outdir, args.bam)
    print(f"Wrote benchmark summary figures to {args.outdir}")


if __name__ == "__main__":
    main()
