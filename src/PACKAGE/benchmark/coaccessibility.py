"""Validate PACKAGE PacBio co-accessibility outputs against a legacy run.

The legacy ranking script emits each constituent-element pair twice (A->B and
B->A), while PACKAGE emits one unordered pair. This module canonicalizes pair IDs
before comparison so agreement statistics describe biological pairs rather than
text-file representation.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class PairRank:
    element_a: str
    element_b: str
    stitched_id: str
    score: float
    rank: int
    is_super: bool

    @property
    def key(self) -> tuple[str, str, str]:
        a, b = sorted((self.element_a, self.element_b))
        return (self.stitched_id, a, b)


@dataclass(frozen=True)
class ClusterRank:
    stitched_id: str
    score: float
    rank: int
    is_super: bool
    counts: tuple[int, int, int, int] | None = None


def _read_pair_ranks(path: Path) -> tuple[dict[tuple[str, str, str], PairRank], int, int]:
    rows: dict[tuple[str, str, str], PairRank] = {}
    n_rows = 0
    duplicate_rows = 0
    with Path(path).open() as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                raise ValueError(f"Expected at least 5 columns at {path}:{line_number}")
            row = PairRank(
                element_a=fields[0],
                element_b=fields[1],
                stitched_id=fields[2],
                score=float(fields[3]),
                rank=int(fields[4]),
                is_super=len(fields) > 5 and fields[5] == "Super",
            )
            n_rows += 1
            existing = rows.get(row.key)
            if existing is not None:
                duplicate_rows += 1
                if not math.isclose(existing.score, row.score, rel_tol=1e-9, abs_tol=1e-12):
                    raise ValueError(
                        f"Directional duplicate has unequal scores for {row.key}: "
                        f"{existing.score} versus {row.score}"
                    )
                continue
            rows[row.key] = row
    return rows, n_rows, duplicate_rows


def _read_cluster_ranks(path: Path) -> dict[str, ClusterRank]:
    rows: dict[str, ClusterRank] = {}
    with Path(path).open() as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                raise ValueError(f"Expected at least 5 columns at {path}:{line_number}")
            row = ClusterRank(
                stitched_id=fields[1],
                score=float(fields[3]),
                rank=int(fields[4]),
                is_super=len(fields) > 5 and fields[5] == "Super",
                counts=(tuple(int(value) for value in fields[6:10]) if len(fields) >= 10 else None),
            )
            if row.stitched_id in rows:
                raise ValueError(f"Duplicate stitched region {row.stitched_id!r} in {path}")
            rows[row.stitched_id] = row
    return rows


def _correlation(x: Iterable[float], y: Iterable[float]) -> float | None:
    x_arr = np.asarray(list(x), dtype=float)
    y_arr = np.asarray(list(y), dtype=float)
    if len(x_arr) < 2 or np.isclose(x_arr.std(), 0) or np.isclose(y_arr.std(), 0):
        return None
    return float(np.corrcoef(x_arr, y_arr)[0, 1])


def _rank_correlation(x: Iterable[float], y: Iterable[float]) -> float | None:
    x_rank = pd.Series(list(x), dtype=float).rank(method="average").to_numpy()
    y_rank = pd.Series(list(y), dtype=float).rank(method="average").to_numpy()
    return _correlation(x_rank, y_rank)


def _classification_counts(
    legacy: Iterable[bool], package: Iterable[bool]
) -> dict[str, int | float | None]:
    legacy_arr = np.asarray(list(legacy), dtype=bool)
    package_arr = np.asarray(list(package), dtype=bool)
    tp = int((legacy_arr & package_arr).sum())
    tn = int((~legacy_arr & ~package_arr).sum())
    fp = int((~legacy_arr & package_arr).sum())
    fn = int((legacy_arr & ~package_arr).sum())
    total = len(legacy_arr)
    return {
        "both_super": tp,
        "both_not_super": tn,
        "package_only_super": fp,
        "legacy_only_super": fn,
        "agreement_fraction": (tp + tn) / total if total else None,
    }


def _top_k_rows(
    legacy: dict[tuple[str, str, str], PairRank],
    package: dict[tuple[str, str, str], PairRank],
    values: list[int],
) -> list[dict[str, int | float]]:
    legacy_order = [row.key for row in sorted(legacy.values(), key=lambda row: row.rank)]
    package_order = [row.key for row in sorted(package.values(), key=lambda row: row.rank)]
    rows: list[dict[str, int | float]] = []
    for requested_k in values:
        k = min(requested_k, len(legacy_order), len(package_order))
        if k == 0:
            continue
        overlap = len(set(legacy_order[:k]) & set(package_order[:k]))
        rows.append(
            {
                "requested_k": requested_k,
                "effective_k": k,
                "overlap": overlap,
                "overlap_fraction": overlap / k,
            }
        )
    return rows


def _fingerprint_lines(path: Path) -> dict[str, Any]:
    """Return an order-independent multiset fingerprint without loading the file."""
    modulus = 1 << 128
    xor_value = 0
    sum_value = 0
    rows = 0
    with Path(path).open("rb") as handle:
        for line in handle:
            digest = int.from_bytes(
                hashlib.blake2b(line.rstrip(b"\n"), digest_size=16).digest(),
                byteorder="big",
            )
            xor_value ^= digest
            sum_value = (sum_value + digest) % modulus
            rows += 1
    return {
        "path": str(path),
        "size_bytes": Path(path).stat().st_size,
        "rows": rows,
        "multiset_xor": f"{xor_value:032x}",
        "multiset_sum": f"{sum_value:032x}",
        "_xor_int": xor_value,
        "_sum_int": sum_value,
    }


def _fingerprint_composes(
    legacy: dict[str, Any], package: dict[str, Any], difference: dict[str, Any]
) -> bool:
    modulus = 1 << 128
    return bool(
        legacy["rows"] == package["rows"] + difference["rows"]
        and legacy["_xor_int"] == (package["_xor_int"] ^ difference["_xor_int"])
        and legacy["_sum_int"]
        == (package["_sum_int"] + difference["_sum_int"]) % modulus
    )


def _public_fingerprint(row: dict[str, Any]) -> dict[str, Any]:
    return {key: value for key, value in row.items() if not key.startswith("_")}


def _object_summary(path: Path, threshold: float) -> dict[str, Any]:
    with Path(path).open() as handle:
        matrix = json.load(handle)
    elements = 0
    cells = 0
    observed = 0
    accessible = 0
    for stitched in matrix:
        for element in stitched.get("enhs", []):
            elements += 1
            values = element.get("fibers", [])
            cells += len(values)
            observed += sum(value is not None for value in values)
            accessible += sum(
                value is not None and float(value) < threshold for value in values
            )
    return {
        "path": str(path),
        "size_bytes": Path(path).stat().st_size,
        "stitched_regions": len(matrix),
        "elements": elements,
        "matrix_cells": cells,
        "observed_cells": observed,
        "accessible_cells": accessible,
    }


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        path.write_text("")
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _write_figures(pair_df: pd.DataFrame, cluster_df: pd.DataFrame, top_k: pd.DataFrame, outdir: Path) -> None:
    try:
        from matplotlib import pyplot as plt
    except ImportError:
        return

    fig, axes = plt.subplots(2, 2, figsize=(10, 8))

    positive = pair_df[(pair_df["legacy_score"] > 0) & (pair_df["package_score"] > 0)]
    axes[0, 0].scatter(
        positive["legacy_score"], positive["package_score"], s=5, alpha=0.25
    )
    if not positive.empty:
        low = min(positive["legacy_score"].min(), positive["package_score"].min())
        high = max(positive["legacy_score"].max(), positive["package_score"].max())
        axes[0, 0].plot([low, high], [low, high], color="black", linewidth=1)
        axes[0, 0].set_xscale("log")
        axes[0, 0].set_yscale("log")
    axes[0, 0].set_xlabel("Legacy corrected score")
    axes[0, 0].set_ylabel("PACKAGE corrected score")
    axes[0, 0].set_title("Constituent-pair scores")

    axes[0, 1].scatter(
        pair_df["legacy_rank"], pair_df["package_rank"], s=5, alpha=0.25
    )
    max_rank = max(pair_df["legacy_rank"].max(), pair_df["package_rank"].max())
    axes[0, 1].plot([1, max_rank], [1, max_rank], color="black", linewidth=1)
    axes[0, 1].set_xscale("log")
    axes[0, 1].set_yscale("log")
    axes[0, 1].set_xlabel("Legacy rank after pair deduplication")
    axes[0, 1].set_ylabel("PACKAGE rank")
    axes[0, 1].set_title("Constituent-pair ranks")

    axes[1, 0].plot(top_k["effective_k"], top_k["overlap_fraction"], marker="o")
    axes[1, 0].set_xscale("log")
    axes[1, 0].set_ylim(0, 1.02)
    axes[1, 0].set_xlabel("Top k pairs")
    axes[1, 0].set_ylabel("Shared fraction")
    axes[1, 0].set_title("Top-k pair overlap")

    positive_clusters = cluster_df[
        (cluster_df["legacy_score"] > 0) & (cluster_df["package_score"] > 0)
    ]
    axes[1, 1].scatter(
        positive_clusters["legacy_score"],
        positive_clusters["package_score"],
        s=5,
        alpha=0.25,
    )
    if not positive_clusters.empty:
        low = min(
            positive_clusters["legacy_score"].min(),
            positive_clusters["package_score"].min(),
        )
        high = max(
            positive_clusters["legacy_score"].max(),
            positive_clusters["package_score"].max(),
        )
        axes[1, 1].plot([low, high], [low, high], color="black", linewidth=1)
        axes[1, 1].set_xscale("log")
        axes[1, 1].set_yscale("log")
    axes[1, 1].set_xlabel("Legacy stitched-region score")
    axes[1, 1].set_ylabel("PACKAGE stitched-region score")
    axes[1, 1].set_title("Shared stitched regions")

    fig.tight_layout()
    fig.savefig(outdir / "coaccess_validation.png", dpi=200)
    fig.savefig(outdir / "coaccess_validation.pdf")
    plt.close(fig)


def run_coaccessibility_validation(
    legacy_ce: Path,
    package_ce: Path,
    legacy_cluster: Path,
    package_cluster: Path,
    outdir: Path,
    *,
    legacy_cov: Path | None = None,
    package_cov: Path | None = None,
    cov_difference: Path | None = None,
    legacy_object: Path | None = None,
    package_object: Path | None = None,
    threshold: float = 0.10,
) -> dict[str, Any]:
    """Compare legacy and PACKAGE co-accessibility results and write figures/tables."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    legacy_pairs, legacy_pair_rows, legacy_duplicates = _read_pair_ranks(legacy_ce)
    package_pairs, package_pair_rows, package_duplicates = _read_pair_ranks(package_ce)
    shared_pair_keys = sorted(set(legacy_pairs) & set(package_pairs))
    pair_rows = []
    for key in shared_pair_keys:
        legacy = legacy_pairs[key]
        package = package_pairs[key]
        pair_rows.append(
            {
                "stitched_id": key[0],
                "element_a": key[1],
                "element_b": key[2],
                "legacy_score": legacy.score,
                "package_score": package.score,
                "legacy_rank": legacy.rank,
                "package_rank": package.rank,
                "legacy_super": legacy.is_super,
                "package_super": package.is_super,
            }
        )
    pair_df = pd.DataFrame(pair_rows)
    pair_df.to_csv(outdir / "pair_rank_agreement.csv", index=False)
    pair_membership_rows = []
    for key in sorted(set(legacy_pairs) | set(package_pairs)):
        legacy = legacy_pairs.get(key)
        package = package_pairs.get(key)
        pair_membership_rows.append(
            {
                "stitched_id": key[0],
                "element_a": key[1],
                "element_b": key[2],
                "status": (
                    "shared" if legacy is not None and package is not None
                    else "legacy_only" if legacy is not None
                    else "package_only"
                ),
                "legacy_score": legacy.score if legacy is not None else "",
                "package_score": package.score if package is not None else "",
                "legacy_rank": legacy.rank if legacy is not None else "",
                "package_rank": package.rank if package is not None else "",
            }
        )
    _write_csv(outdir / "pair_membership.csv", pair_membership_rows)

    legacy_clusters = _read_cluster_ranks(legacy_cluster)
    package_clusters = _read_cluster_ranks(package_cluster)
    shared_cluster_ids = sorted(set(legacy_clusters) & set(package_clusters))
    cluster_rows = []
    for stitched_id in shared_cluster_ids:
        legacy = legacy_clusters[stitched_id]
        package = package_clusters[stitched_id]
        cluster_rows.append(
            {
                "stitched_id": stitched_id,
                "legacy_score": legacy.score,
                "package_score": package.score,
                "legacy_rank": legacy.rank,
                "package_rank": package.rank,
                "legacy_super": legacy.is_super,
                "package_super": package.is_super,
            }
        )
    cluster_df = pd.DataFrame(cluster_rows)
    cluster_df.to_csv(outdir / "cluster_rank_agreement.csv", index=False)
    cluster_membership_rows = []
    package_only_half_jaccard = 0
    for stitched_id in sorted(set(legacy_clusters) | set(package_clusters)):
        legacy = legacy_clusters.get(stitched_id)
        package = package_clusters.get(stitched_id)
        status = (
            "shared" if legacy is not None and package is not None
            else "legacy_only" if legacy is not None
            else "package_only"
        )
        counts = package.counts if package is not None else None
        half_jaccard_candidate = bool(
            status == "package_only"
            and counts is not None
            and 2 * counts[0] == counts[1] + counts[2]
        )
        package_only_half_jaccard += int(half_jaccard_candidate)
        cluster_membership_rows.append(
            {
                "stitched_id": stitched_id,
                "status": status,
                "legacy_score": legacy.score if legacy is not None else "",
                "package_score": package.score if package is not None else "",
                "legacy_rank": legacy.rank if legacy is not None else "",
                "package_rank": package.rank if package is not None else "",
                "package_best_pair_a": counts[0] if counts is not None else "",
                "package_best_pair_b": counts[1] if counts is not None else "",
                "package_best_pair_c": counts[2] if counts is not None else "",
                "package_best_pair_d": counts[3] if counts is not None else "",
                "legacy_half_jaccard_sentinel_candidate": half_jaccard_candidate,
            }
        )
    _write_csv(outdir / "cluster_membership.csv", cluster_membership_rows)

    top_k_rows = _top_k_rows(
        legacy_pairs,
        package_pairs,
        [10, 50, 100, 500, 1000, 5000, 10_000, 50_000],
    )
    _write_csv(outdir / "top_k_overlap.csv", top_k_rows)
    top_k_df = pd.DataFrame(top_k_rows)

    pair_summary = {
        "legacy_rows": legacy_pair_rows,
        "legacy_unique_unordered_pairs": len(legacy_pairs),
        "legacy_directional_duplicate_rows": legacy_duplicates,
        "package_rows": package_pair_rows,
        "package_unique_unordered_pairs": len(package_pairs),
        "package_duplicate_rows": package_duplicates,
        "shared_pairs": len(shared_pair_keys),
        "legacy_only_pairs": len(set(legacy_pairs) - set(package_pairs)),
        "package_only_pairs": len(set(package_pairs) - set(legacy_pairs)),
        "score_pearson": _correlation(
            pair_df["legacy_score"], pair_df["package_score"]
        ),
        "score_spearman": _rank_correlation(
            pair_df["legacy_score"], pair_df["package_score"]
        ),
        "rank_spearman": _rank_correlation(
            pair_df["legacy_rank"], pair_df["package_rank"]
        ),
        "super_classification": _classification_counts(
            pair_df["legacy_super"], pair_df["package_super"]
        ),
    }
    cluster_summary = {
        "legacy_clusters": len(legacy_clusters),
        "package_clusters": len(package_clusters),
        "shared_clusters": len(shared_cluster_ids),
        "legacy_only_clusters": len(set(legacy_clusters) - set(package_clusters)),
        "package_only_clusters": len(set(package_clusters) - set(legacy_clusters)),
        "package_only_half_jaccard_sentinel_candidates": package_only_half_jaccard,
        "score_pearson": _correlation(
            cluster_df["legacy_score"], cluster_df["package_score"]
        ),
        "score_spearman": _rank_correlation(
            cluster_df["legacy_score"], cluster_df["package_score"]
        ),
        "rank_spearman": _rank_correlation(
            cluster_df["legacy_rank"], cluster_df["package_rank"]
        ),
        "super_classification": _classification_counts(
            cluster_df["legacy_super"], cluster_df["package_super"]
        ),
    }

    summary: dict[str, Any] = {
        "threshold": threshold,
        "comparison_notes": {
            "legacy_pair_rows_are_directional": True,
            "package_pair_rows_are_unordered": True,
            "legacy_cluster_code_uses_0.5_as_invalid_pair_sentinel": True,
            "legacy_distance_fit_uses_directional_duplicates": True,
        },
        "pairs": pair_summary,
        "clusters": cluster_summary,
        "top_k_overlap": top_k_rows,
    }

    if legacy_cov is not None and package_cov is not None:
        legacy_fp = _fingerprint_lines(legacy_cov)
        package_fp = _fingerprint_lines(package_cov)
        cov_summary: dict[str, Any] = {
            "legacy": _public_fingerprint(legacy_fp),
            "package": _public_fingerprint(package_fp),
            "exact_multiset_match": (
                legacy_fp["rows"] == package_fp["rows"]
                and legacy_fp["_xor_int"] == package_fp["_xor_int"]
                and legacy_fp["_sum_int"] == package_fp["_sum_int"]
            ),
        }
        if cov_difference is not None:
            difference_fp = _fingerprint_lines(cov_difference)
            cov_summary["known_difference"] = _public_fingerprint(difference_fp)
            cov_summary["package_plus_difference_matches_legacy"] = _fingerprint_composes(
                legacy_fp, package_fp, difference_fp
            )
        summary["cov"] = cov_summary

    if legacy_object is not None and package_object is not None:
        summary["objects"] = {
            "legacy": _object_summary(legacy_object, threshold),
            "package": _object_summary(package_object, threshold),
        }

    (outdir / "validation_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    _write_figures(pair_df, cluster_df, top_k_df, outdir)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--legacy-ce", type=Path, required=True)
    parser.add_argument("--package-ce", type=Path, required=True)
    parser.add_argument("--legacy-cluster", type=Path, required=True)
    parser.add_argument("--package-cluster", type=Path, required=True)
    parser.add_argument("--legacy-cov", type=Path)
    parser.add_argument("--package-cov", type=Path)
    parser.add_argument("--cov-difference", type=Path)
    parser.add_argument("--legacy-object", type=Path)
    parser.add_argument("--package-object", type=Path)
    parser.add_argument("--threshold", type=float, default=0.10)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()

    summary = run_coaccessibility_validation(
        args.legacy_ce,
        args.package_ce,
        args.legacy_cluster,
        args.package_cluster,
        args.outdir,
        legacy_cov=args.legacy_cov,
        package_cov=args.package_cov,
        cov_difference=args.cov_difference,
        legacy_object=args.legacy_object,
        package_object=args.package_object,
        threshold=args.threshold,
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
