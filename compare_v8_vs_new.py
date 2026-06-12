#!/usr/bin/env python3
"""Validate the new PACKAGE.db.database port against the V8 reference implementation.

Compares outputs of the V8 ``FiberDatabase`` and the new ``PACKAGE.db.FiberDatabase``
on the same HDF5 database file, across multiple query types. Both implementations
should produce identical results (modulo row ordering and float precision).

Usage
=====

    # On HPC, with both V8 code and the new package available:
    python compare_v8_vs_new.py \\
        --v8-code-dir /nfs/roberts/scratch/pi_zs296/mc3365/long_read/first_run/code/build_db \\
        --db /nfs/roberts/scratch/pi_zs296/mc3365/long_read/first_run/out/h5/fiber_database_v8.h5 \\
        --sample d0 \\
        --annotation CGI \\
        --max-regions 50

Exit code is 0 if all checks pass, 1 otherwise. Diffs are printed to stderr.
"""

from __future__ import annotations

import argparse
import sys
import time
import traceback
from pathlib import Path

import pandas as pd


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

GREEN = "\033[92m"
RED = "\033[91m"
YELLOW = "\033[93m"
RESET = "\033[0m"


def ok(msg: str) -> None:
    print(f"{GREEN}✓{RESET} {msg}")


def fail(msg: str) -> None:
    print(f"{RED}✗{RESET} {msg}", file=sys.stderr)


def warn(msg: str) -> None:
    print(f"{YELLOW}!{RESET} {msg}")


def info(msg: str) -> None:
    print(f"  {msg}")


def sort_df_for_compare(df: pd.DataFrame) -> pd.DataFrame:
    """Canonicalize a result DataFrame so order-of-rows differences don't matter.

    Sorts by all columns that look like row identifiers; resets the index so
    pd.testing.assert_frame_equal will compare on values, not on index labels.
    """
    sort_cols = [
        c for c in ("region_id", "annotation", "sample", "fiber_id") if c in df.columns
    ]
    if not sort_cols:
        return df.reset_index(drop=True)
    return df.sort_values(sort_cols).reset_index(drop=True)


def compare_dataframes(
    df_v8: pd.DataFrame, df_new: pd.DataFrame, label: str
) -> bool:
    """Compare two DataFrames after canonicalizing row order.

    Returns True if equal (allowing small float differences), False otherwise.
    Prints a diagnostic dump on failure.
    """
    df_v8_s = sort_df_for_compare(df_v8)
    df_new_s = sort_df_for_compare(df_new)

    # Shape check first; gives most useful error if column sets differ.
    if df_v8_s.shape != df_new_s.shape:
        fail(f"{label}: shape mismatch — V8={df_v8_s.shape} vs new={df_new_s.shape}")
        info(f"V8 columns: {sorted(df_v8_s.columns)}")
        info(f"New columns: {sorted(df_new_s.columns)}")
        # Show first rows of each so you can eyeball the difference
        info(f"V8 head:\n{df_v8_s.head(3)}")
        info(f"New head:\n{df_new_s.head(3)}")
        return False

    # Column set check
    if set(df_v8_s.columns) != set(df_new_s.columns):
        fail(f"{label}: column set mismatch")
        info(f"V8 only: {set(df_v8_s.columns) - set(df_new_s.columns)}")
        info(f"New only: {set(df_new_s.columns) - set(df_v8_s.columns)}")
        return False

    # Reorder new columns to match V8 column order so assert_frame_equal compares
    # by value, not column position
    df_new_s = df_new_s[df_v8_s.columns]

    try:
        pd.testing.assert_frame_equal(
            df_v8_s,
            df_new_s,
            check_dtype=False,   # int64 vs Int64 etc. are equivalent for our purposes
            check_exact=False,   # allow tiny float precision differences
            rtol=1e-6,
            atol=1e-9,
        )
        ok(f"{label}: {len(df_v8_s)} rows identical")
        return True
    except AssertionError as e:
        fail(f"{label}: DataFrames differ")
        info(str(e).split("\n")[0])  # one-line summary
        # Show first few differing rows for debugging
        try:
            diff_mask = (df_v8_s != df_new_s).any(axis=1)
            if diff_mask.any():
                info(f"First {min(5, diff_mask.sum())} differing rows:")
                info(f"V8:\n{df_v8_s[diff_mask].head()}")
                info(f"New:\n{df_new_s[diff_mask].head()}")
        except Exception:
            pass  # diff logic itself can fail; don't mask the real assertion error
        return False


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("Usage")[0])
    parser.add_argument(
        "--v8-code-dir",
        type=Path,
        required=True,
        help="Directory containing fiber_database_v8.py",
    )
    parser.add_argument(
        "--db", type=Path, required=True, help="Path to an existing V8 HDF5 database."
    )
    parser.add_argument(
        "--sample",
        default=None,
        help="Sample name to query (defaults to first sample in db).",
    )
    parser.add_argument(
        "--annotation",
        default="CGI",
        help="Annotation type to query (default: CGI).",
    )
    parser.add_argument(
        "--max-regions",
        type=int,
        default=50,
        help="Limit each query to this many regions for speed (default: 50).",
    )
    args = parser.parse_args()

    if not args.v8_code_dir.exists():
        fail(f"V8 code directory not found: {args.v8_code_dir}")
        return 1
    if not args.db.exists():
        fail(f"Database not found: {args.db}")
        return 1

    # ----- Import both implementations -----
    sys.path.insert(0, str(args.v8_code_dir))
    try:
        from fiber_database_v8 import FiberDatabase as V8FiberDatabase
    except ImportError:
        fail(f"Could not import fiber_database_v8 from {args.v8_code_dir}")
        traceback.print_exc()
        return 1
    try:
        from PACKAGE.db import FiberDatabase as NewFiberDatabase
    except ImportError:
        fail("Could not import PACKAGE.db.FiberDatabase — is the new package installed?")
        traceback.print_exc()
        return 1

    ok("Imports succeeded for both V8 and new")

    # ----- Open both -----
    print(f"\nOpening database: {args.db}")
    db_v8 = V8FiberDatabase(str(args.db))
    db_new = NewFiberDatabase(args.db)
    sample = args.sample or db_v8.samples[0]
    info(f"Using sample: {sample}")
    info(f"Samples found V8: {db_v8.samples}, new: {db_new.samples}")
    if db_v8.samples != db_new.samples:
        fail("Sample lists differ")
        return 1
    ok("Sample lists match")

    failures = 0

    # ----- Test 1: list_annotations -----
    print("\n[1/6] list_annotations()")
    a_v8 = db_v8.list_annotations()
    a_new = db_new.list_annotations()
    if a_v8 == a_new:
        ok(f"Annotation lists identical ({len(a_v8)} types)")
    else:
        fail("Annotation lists differ")
        info(f"V8 only: {set(a_v8) - set(a_new)}")
        info(f"New only: {set(a_new) - set(a_v8)}")
        failures += 1

    # ----- Test 2: get_chromosomes -----
    print("\n[2/6] get_chromosomes()")
    c_v8 = db_v8.get_chromosomes(sample)
    c_new = db_new.get_chromosomes(sample)
    if c_v8 == c_new:
        ok(f"Chromosome lists identical ({len(c_v8)} chroms)")
    else:
        fail("Chromosome lists differ")
        info(f"V8: {c_v8}")
        info(f"New: {c_new}")
        failures += 1

    # ----- Test 3: get_summary -----
    print("\n[3/6] get_summary()")
    s_v8 = db_v8.get_summary()
    s_new = db_new.get_summary()
    # Both should have version, genome, samples, n_annotations, plus n_fibers_X per sample.
    # Compare just the shared keys.
    shared = set(s_v8) & set(s_new)
    diffs = {k: (s_v8[k], s_new[k]) for k in shared if s_v8[k] != s_new[k]}
    if not diffs:
        ok(f"Summaries match on {len(shared)} shared keys")
    else:
        fail("Summary values differ")
        for k, (v8, new) in diffs.items():
            info(f"  {k}: V8={v8}, new={new}")
        failures += 1

    # ----- Test 4: bulk query — vectorized query_annotation_fast -----
    print(f"\n[4/6] query_annotation_fast('{args.annotation}', max_regions={args.max_regions})")
    t0 = time.time()
    df_v8 = db_v8.query_annotation_fast(
        args.annotation, sample=sample, max_regions=args.max_regions
    )
    t_v8 = time.time() - t0
    t0 = time.time()
    df_new = db_new.query_annotation_fast(
        args.annotation, sample=sample, max_regions=args.max_regions
    )
    t_new = time.time() - t0
    info(f"V8 time: {t_v8:.2f}s, new time: {t_new:.2f}s")
    if not compare_dataframes(df_v8, df_new, "query_annotation_fast (default layers)"):
        failures += 1

    # ----- Test 5: bulk query with all 5 layers -----
    print(f"\n[5/6] query_annotation_fast with all layers")
    all_layers = ["nucleosomes", "5mC", "5hmC", "6mA", "msp"]
    df_v8_all = db_v8.query_annotation_fast(
        args.annotation,
        sample=sample,
        feature_types=all_layers,
        max_regions=args.max_regions,
    )
    df_new_all = db_new.query_annotation_fast(
        args.annotation,
        sample=sample,
        feature_types=all_layers,
        max_regions=args.max_regions,
    )
    if not compare_dataframes(df_v8_all, df_new_all, "query_annotation_fast (all layers)"):
        failures += 1

    # ----- Test 6: per-fiber accessors on one fiber -----
    print(f"\n[6/6] per-fiber accessors on one fiber")
    # Pick a real fiber from the first chromosome that has data
    chrom0 = c_v8[0]
    fids_v8 = db_v8.get_fibers_at(chrom0, 0, 10_000_000, sample=sample)
    fids_new = db_new.get_fibers_at(chrom0, 0, 10_000_000, sample=sample)
    # Note: order may differ between spatial-index path and array scan; compare as sets
    if set(fids_v8) == set(fids_new):
        ok(f"get_fibers_at: {len(fids_v8)} fibers, identical sets")
    else:
        only_v8 = set(fids_v8) - set(fids_new)
        only_new = set(fids_new) - set(fids_v8)
        fail(f"get_fibers_at: {len(only_v8)} V8-only, {len(only_new)} new-only")
        if only_v8:
            info(f"V8 only (first 3): {list(only_v8)[:3]}")
        if only_new:
            info(f"New only (first 3): {list(only_new)[:3]}")
        failures += 1

    if fids_v8:
        test_fiber = fids_v8[0]
        for method_name in ("get_nucleosomes", "get_methylation", "get_msp"):
            v8_method = getattr(db_v8, method_name)
            new_method = getattr(db_new, method_name)
            try:
                if method_name == "get_methylation":
                    res_v8 = v8_method(test_fiber, chrom0, "5mC", sample)
                    res_new = new_method(test_fiber, chrom0, "5mC", sample)
                else:
                    res_v8 = v8_method(test_fiber, chrom0, sample)
                    res_new = new_method(test_fiber, chrom0, sample)
                # Both should be dicts with same keys; compare each numpy array.
                if set(res_v8.keys()) != set(res_new.keys()):
                    fail(f"{method_name}: key set differs")
                    info(f"V8 keys: {sorted(res_v8.keys())}")
                    info(f"New keys: {sorted(res_new.keys())}")
                    failures += 1
                    continue
                all_match = all(
                    (res_v8[k] == res_new[k]).all() if hasattr(res_v8[k], "all")
                    else res_v8[k] == res_new[k]
                    for k in res_v8
                )
                if all_match:
                    ok(f"{method_name}: outputs identical")
                else:
                    fail(f"{method_name}: values differ")
                    failures += 1
            except Exception as e:
                fail(f"{method_name}: exception raised — {e}")
                failures += 1

    db_v8.close()
    db_new.close()

    # ----- Summary -----
    print()
    if failures == 0:
        ok("ALL CHECKS PASSED — V8 and new produce identical outputs.")
        return 0
    else:
        fail(f"{failures} checks failed — see above.")
        return 1


if __name__ == "__main__":
    sys.exit(main())
