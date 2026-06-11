"""Database builder: extracted intermediate files → HDF5 database.

Entry point called by ``PACKAGE build``. Reads the intermediate files produced by
``PACKAGE extract`` (or by manual modkit/fibertools-rs invocations) and writes them
into the HDF5 schema defined in :mod:`PACKAGE.db.schema`.

This is a refactor of ``build_fiber_database_v8.py`` with the following changes from V8:

  - File paths come from a :class:`PACKAGE.config.Config` rather than being derived
    from a ``base_dir`` and sample-prefix convention. Each layer's file path is
    explicit in the YAML, giving full flexibility.
  - The set of layers built per sample is configurable. Only the ``nucleosomes``
    layer is required (it defines which fibers exist on each chromosome);
    ``5mC``, ``5hmC``, ``6mA``, and ``msp`` are optional and skipped if their
    paths aren't listed in the sample config.
  - All HDF5 paths come from :mod:`PACKAGE.db.schema` instead of inline string
    literals (matches the convention established by the query layer).
  - Logging via :func:`PACKAGE.utils.get_logger` replaces V8's print statements.

The builder is split across three batches matching the V8 source structure:

  - Batch 1 (this file currently): foundation + format parsers
  - Batch 2: per-chromosome assembly helpers (ID mapping, sorting, writing)
  - Batch 3: orchestration (per-sample build + annotation loading + main entry)
"""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from tqdm import tqdm

from PACKAGE.config import Config, SampleConfig
from PACKAGE.db import schema
from PACKAGE.utils import get_logger, smart_open

log = get_logger(__name__)


# ---------------------------------------------------------------------------
# Public entry point (called by the CLI)
# ---------------------------------------------------------------------------

def build_database(config: Config, samples: list[str] | None = None) -> None:
    """Build an HDF5 database from extracted intermediate files.

    Args:
        config: Validated PACKAGE config.
        samples: Subset of sample names. If None, builds all in config.
    """
    samples = samples or config.sample_names
    log.info(f"Build target: {config.output_path}")
    log.info(f"Samples to build: {samples}")

    # Resolve each requested sample name to its SampleConfig (raises KeyError if missing).
    sample_configs = [config.get_sample(name) for name in samples]

    builder = FiberDatabaseBuilder(config, sample_configs)
    builder.build()


# ---------------------------------------------------------------------------
# Builder class
# ---------------------------------------------------------------------------

class FiberDatabaseBuilder:
    """Internal builder. Most users should call :func:`build_database` instead.

    The class is exposed mainly so tests can swap in mocked file readers, or so
    advanced users can subclass to customize specific steps.
    """

    # Required column count for modkit extract TSV. Used as a sanity-check filter.
    _MIN_MODKIT_COLS: int = 14

    # Column indices in modkit's `extract full` output (1-based in V8 comments,
    # 0-based here as Python expects). These are fixed by modkit's spec, so we
    # name them as constants instead of magic numbers in the parser.
    _MODKIT_COL_READ_ID = 0
    _MODKIT_COL_REF_POS = 2
    _MODKIT_COL_CHROM = 3
    _MODKIT_COL_MOD_QUAL = 12
    _MODKIT_COL_MOD_CODE = 13

    def __init__(
        self,
        config: Config,
        sample_configs: list[SampleConfig],
    ) -> None:
        self.config = config
        self.sample_configs = sample_configs
        self.output_path = config.output_path
        self.methylation_threshold = config.parameters.methylation_threshold
        self.min_msp_size = config.parameters.min_msp_size

    # ------------------------------------------------------------------
    # Input validation
    # ------------------------------------------------------------------

    def validate_inputs(self) -> None:
        """Check that all input files referenced in the config exist.

        Catches typos and stale paths at job start rather than after hours of
        processing. Raises FileNotFoundError listing every missing file.
        """
        missing: list[str] = []
        for sc in self.sample_configs:
            for layer, path in sc.layers.items():
                if not path.exists():
                    missing.append(f"{sc.name}.{layer}: {path}")
        # Annotation files too.
        for name, path in self.config.annotations.items():
            if not path.exists():
                missing.append(f"annotation '{name}': {path}")
        if missing:
            raise FileNotFoundError(
                "The following input files do not exist:\n  - "
                + "\n  - ".join(missing)
            )

    # ------------------------------------------------------------------
    # Format parsers
    # ------------------------------------------------------------------
    # These are platform-agnostic parsers for the intermediate file formats
    # produced by modkit + fibertools-rs ft extract. Each yields per-call records
    # filtered by a per-chromosome "known fibers" set so we don't carry around
    # data for reads we'll later discard.

    def parse_bed12_single_base(
        self,
        filepath: Path,
        known_fibers: dict[str, set[bytes]],
    ) -> dict[str, dict[str, list]]:
        """Parse a BED12 file where each block is one single-base call.

        This is the format produced by ``ft extract --m6a`` (one block per called
        6mA position, plus first/last blocks as boundary markers which we skip).
        Each block has size 1; position = ``abs_start + block_offset``.

        Args:
            filepath: BED12 file path (may be .gz; smart_open handles that).
            known_fibers: ``{chrom: {fiber_id_bytes}}``. Reads not in their
                chromosome's set are dropped.

        Returns:
            ``{chrom: {"fids": [bytes,...], "pos": [int,...]}}``
        """
        result: dict[str, dict[str, list]] = defaultdict(
            lambda: {"fids": [], "pos": []}
        )
        with smart_open(filepath) as f:
            for line in tqdm(f, desc=f"  {filepath.name}", unit_scale=True):
                cols = line.rstrip("\n").split("\t")
                if len(cols) < 12:
                    continue
                chrom = cols[0]
                if chrom not in known_fibers:
                    continue
                read_id = cols[3].encode()
                if read_id not in known_fibers[chrom]:
                    continue
                abs_start = int(cols[1])
                offsets = [int(x) for x in cols[11].rstrip(",").split(",")]
                # Skip first and last blocks (ft extract boundary markers).
                for offset in offsets[1:-1]:
                    result[chrom]["fids"].append(read_id)
                    result[chrom]["pos"].append(abs_start + offset)
        return dict(result)

    def parse_bed12_intervals(
        self,
        filepath: Path,
        known_fibers: dict[str, set[bytes]],
        min_size: int = 0,
    ) -> dict[str, dict[str, list]]:
        """Parse a BED12 file where each block is an interval (start, end).

        This is the format produced by ``ft extract --msp``. Each block has a
        nonzero size; we filter blocks smaller than ``min_size`` as artifacts.
        First and last blocks (boundary markers) are skipped.

        Returns:
            ``{chrom: {"fids": [...], "starts": [...], "ends": [...], "widths": [...]}}``
        """
        result: dict[str, dict[str, list]] = defaultdict(
            lambda: {"fids": [], "starts": [], "ends": [], "widths": []}
        )
        with smart_open(filepath) as f:
            for line in tqdm(f, desc=f"  {filepath.name}", unit_scale=True):
                cols = line.rstrip("\n").split("\t")
                if len(cols) < 12:
                    continue
                chrom = cols[0]
                if chrom not in known_fibers:
                    continue
                read_id = cols[3].encode()
                if read_id not in known_fibers[chrom]:
                    continue
                abs_start = int(cols[1])
                sizes = [int(x) for x in cols[10].rstrip(",").split(",")]
                offsets = [int(x) for x in cols[11].rstrip(",").split(",")]
                # Skip first and last (boundary markers).
                for size, offset in zip(sizes[1:-1], offsets[1:-1], strict=False):
                    if size < min_size:
                        continue
                    s = abs_start + offset
                    e = s + size
                    result[chrom]["fids"].append(read_id)
                    result[chrom]["starts"].append(s)
                    result[chrom]["ends"].append(e)
                    result[chrom]["widths"].append(size)
        return dict(result)

    def parse_modkit_extract(
        self,
        filepath: Path,
        known_fibers: dict[str, set[bytes]],
        want_5mc: bool = True,
        want_5hmc: bool = True,
    ) -> tuple[dict[str, dict[str, list]], dict[str, dict[str, list]]]:
        """Parse modkit extract TSV; route rows to 5mC and/or 5hmC buckets.

        modkit extract produces one row per (read, CpG, mod_type) combination. The
        same reference position appears twice: once with mod_code='m' (5mC) and
        once with mod_code='h' (5hmC). This parser reads the file once and routes
        each row to the appropriate output based on ``mod_code``.

        Args:
            filepath: modkit extract TSV (may be .gz).
            known_fibers: per-chromosome fiber filter (same shape as other parsers).
            want_5mc, want_5hmc: skip rows whose mod_code corresponds to a layer
                the caller doesn't want. This avoids parsing-time work when a user
                opts out of one of the methylation channels.

        Returns:
            Tuple ``(mc_data, hmc_data)``, each a dict with the same structure
            as the BED12 parsers: ``{chrom: {"fids": [...], "pos": [...], "prob": [...]}}``.
            If a layer isn't wanted, its dict is empty.
        """
        mc_data: dict[str, dict[str, list]] = defaultdict(
            lambda: {"fids": [], "pos": [], "prob": []}
        )
        hmc_data: dict[str, dict[str, list]] = defaultdict(
            lambda: {"fids": [], "pos": [], "prob": []}
        )

        with smart_open(filepath) as f:
            # Read and discard the header.
            _header = f.readline()
            for line in tqdm(f, desc=f"  {filepath.name}", unit_scale=True):
                cols = line.split("\t")
                if len(cols) < self._MIN_MODKIT_COLS:
                    continue
                try:
                    read_id = cols[self._MODKIT_COL_READ_ID]
                    ref_pos = int(cols[self._MODKIT_COL_REF_POS])
                    chrom = cols[self._MODKIT_COL_CHROM]
                    mod_qual = float(cols[self._MODKIT_COL_MOD_QUAL])
                    mod_code = cols[self._MODKIT_COL_MOD_CODE].rstrip()
                except (ValueError, IndexError):
                    # Skip malformed rows silently (matches V8 behavior).
                    continue
                if chrom not in known_fibers:
                    continue
                rid_bytes = read_id.encode()
                if rid_bytes not in known_fibers[chrom]:
                    continue
                if mod_code == "m" and want_5mc:
                    mc_data[chrom]["fids"].append(rid_bytes)
                    mc_data[chrom]["pos"].append(ref_pos)
                    mc_data[chrom]["prob"].append(mod_qual)
                elif mod_code == "h" and want_5hmc:
                    hmc_data[chrom]["fids"].append(rid_bytes)
                    hmc_data[chrom]["pos"].append(ref_pos)
                    hmc_data[chrom]["prob"].append(mod_qual)
        return dict(mc_data), dict(hmc_data)

    # ------------------------------------------------------------------
    # Integer ID mapping and per-chromosome assembly
    # ------------------------------------------------------------------

    @staticmethod
    def _encode_fiber_id(fiber_id: str | bytes | np.bytes_) -> bytes:
        """Normalize fiber IDs to bytes for HDF5 string tables and dict keys."""
        if isinstance(fiber_id, bytes):
            return fiber_id
        if isinstance(fiber_id, np.bytes_):
            return bytes(fiber_id)
        return str(fiber_id).encode()

    def _build_fiber_int_map(
        self,
        fiber_ids_sorted: list[str | bytes] | np.ndarray,
    ) -> tuple[dict[bytes, np.uint32], np.ndarray]:
        """Build bytes UUID -> uint32 mapping for one chromosome.

        This ports V8's per-chromosome integer ID scheme. The returned table is
        stored as ``fiber_id_table``; row index N is the original fiber ID for
        integer ID N.
        """
        table = np.array(
            [self._encode_fiber_id(fid) for fid in fiber_ids_sorted],
            dtype="S50",
        )
        str_to_int = {fid: np.uint32(i) for i, fid in enumerate(table)}
        return str_to_int, table

    def _map_fids_to_int(
        self,
        fid_strings: list[str | bytes] | np.ndarray | pd.Series,
        str_to_int: dict[bytes, np.uint32],
    ) -> np.ndarray:
        """Convert fiber IDs to uint32 using a per-chromosome mapping."""
        return np.array(
            [str_to_int[self._encode_fiber_id(fid)] for fid in fid_strings],
            dtype=np.uint32,
        )

    def _sort_and_assemble(
        self,
        fids_list: list[str | bytes] | np.ndarray,
        arrays_dict: dict[str, list | np.ndarray],
        str_to_int: dict[bytes, np.uint32],
    ) -> tuple[np.ndarray, dict[str, np.ndarray]]:
        """Sort parallel arrays by uint32 fiber ID.

        Stable sorting preserves within-fiber input order, matching V8.
        """
        if len(fids_list) == 0:
            empty_fids = np.array([], dtype=np.uint32)
            return empty_fids, {
                name: np.array([], dtype=np.int32)
                for name in arrays_dict
            }

        int_fids = self._map_fids_to_int(fids_list, str_to_int)
        order = np.argsort(int_fids, kind="stable")
        sorted_fids = int_fids[order]
        sorted_arrays = {}
        for name, values in arrays_dict.items():
            arr = np.array(values)
            sorted_arrays[name] = arr[order]
        return sorted_fids, sorted_arrays

    def load_nucleosome_inputs(
        self,
        sample_config: SampleConfig,
    ) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, set[bytes]]]:
        """Load nucleosome CSV and derive fiber metadata / known fiber sets.

        Nucleosomes are required because they define the set of fibers stored in
        the database, exactly as in V8.
        """
        nuc_path = sample_config.layers["nucleosomes"]
        nuc_df = pd.read_csv(nuc_path)
        nuc_df["linker_dist"] = nuc_df["linker_dist"].fillna(-1)
        fiber_meta = (
            nuc_df.groupby(["chrom", "read_id"])
            .agg(start=("nuc_start", "min"), end=("nuc_end", "max"))
            .reset_index()
            .sort_values(["chrom", "read_id"])
            .reset_index(drop=True)
        )
        known_fibers: dict[str, set[bytes]] = {}
        for chrom in sorted(nuc_df["chrom"].unique()):
            ids = fiber_meta.loc[fiber_meta["chrom"] == chrom, "read_id"]
            known_fibers[chrom] = {self._encode_fiber_id(fid) for fid in ids}
        return nuc_df, fiber_meta, known_fibers

    def assemble_chromosome_data(
        self,
        chrom: str,
        fiber_meta: pd.DataFrame,
        nuc_df: pd.DataFrame,
        mc_data: dict[str, dict[str, list]] | None = None,
        hmc_data: dict[str, dict[str, list]] | None = None,
        ma_data: dict[str, dict[str, list]] | None = None,
        msp_data: dict[str, dict[str, list]] | None = None,
    ) -> dict[str, np.ndarray | int]:
        """Assemble all arrays for one chromosome in schema-compatible order."""
        mc_data = mc_data or {}
        hmc_data = hmc_data or {}
        ma_data = ma_data or {}
        msp_data = msp_data or {}

        fm = fiber_meta[fiber_meta["chrom"] == chrom].sort_values("read_id")
        nuc_chrom = nuc_df[nuc_df["chrom"] == chrom]

        unique_fids_sorted = sorted(fm["read_id"].unique())
        str_to_int, fiber_id_table = self._build_fiber_int_map(unique_fids_sorted)

        fm_int_ids = self._map_fids_to_int(fm["read_id"].values, str_to_int)

        nuc_int_fids = self._map_fids_to_int(nuc_chrom["read_id"].values, str_to_int)
        nuc_order = np.argsort(nuc_int_fids, kind="stable")

        mc_int_fids, mc_sorted = self._sort_and_assemble(
            mc_data.get(chrom, {}).get("fids", []),
            {
                "pos": mc_data.get(chrom, {}).get("pos", []),
                "prob": mc_data.get(chrom, {}).get("prob", []),
            },
            str_to_int,
        )
        mc_prob = (
            mc_sorted.get("prob", np.array([], dtype=np.float32)).astype(np.float32)
            if len(mc_int_fids) > 0
            else np.array([], dtype=np.float32)
        )
        mc_is_meth = (
            (mc_prob >= self.methylation_threshold).astype(np.bool_)
            if len(mc_prob) > 0
            else np.array([], dtype=np.bool_)
        )

        hmc_int_fids, hmc_sorted = self._sort_and_assemble(
            hmc_data.get(chrom, {}).get("fids", []),
            {
                "pos": hmc_data.get(chrom, {}).get("pos", []),
                "prob": hmc_data.get(chrom, {}).get("prob", []),
            },
            str_to_int,
        )
        hmc_prob = (
            hmc_sorted.get("prob", np.array([], dtype=np.float32)).astype(np.float32)
            if len(hmc_int_fids) > 0
            else np.array([], dtype=np.float32)
        )
        hmc_is_meth = (
            (hmc_prob >= self.methylation_threshold).astype(np.bool_)
            if len(hmc_prob) > 0
            else np.array([], dtype=np.bool_)
        )

        ma_int_fids, ma_sorted = self._sort_and_assemble(
            ma_data.get(chrom, {}).get("fids", []),
            {"pos": ma_data.get(chrom, {}).get("pos", [])},
            str_to_int,
        )

        msp_int_fids, msp_sorted = self._sort_and_assemble(
            msp_data.get(chrom, {}).get("fids", []),
            {
                "starts": msp_data.get(chrom, {}).get("starts", []),
                "ends": msp_data.get(chrom, {}).get("ends", []),
                "widths": msp_data.get(chrom, {}).get("widths", []),
            },
            str_to_int,
        )

        return {
            "n_fibers": len(fm),
            "fiber_id_table": fiber_id_table,
            "fm_int_ids": fm_int_ids,
            "fib_starts": fm["start"].values.astype(np.uint32),
            "fib_ends": fm["end"].values.astype(np.uint32),
            "nuc_int_ids": nuc_int_fids[nuc_order],
            "nuc_starts": nuc_chrom["nuc_start"].values[nuc_order].astype(np.uint32),
            "nuc_ends": nuc_chrom["nuc_end"].values[nuc_order].astype(np.uint32),
            "nuc_widths": nuc_chrom["width"].values[nuc_order].astype(np.uint32),
            "nuc_linkers": nuc_chrom["linker_dist"].values[nuc_order].astype(np.float32),
            "mc_int_ids": mc_int_fids,
            "mc_pos": mc_sorted.get("pos", np.array([], dtype=np.int32)).astype(np.uint32),
            "mc_prob": mc_prob,
            "mc_is_meth": mc_is_meth,
            "hmc_int_ids": hmc_int_fids,
            "hmc_pos": hmc_sorted.get("pos", np.array([], dtype=np.int32)).astype(np.uint32),
            "hmc_prob": hmc_prob,
            "hmc_is_meth": hmc_is_meth,
            "ma_int_ids": ma_int_fids,
            "ma_pos": ma_sorted.get("pos", np.array([], dtype=np.int32)).astype(np.uint32),
            "msp_int_ids": msp_int_fids,
            "msp_starts": msp_sorted.get("starts", np.array([], dtype=np.int32)).astype(np.uint32),
            "msp_ends": msp_sorted.get("ends", np.array([], dtype=np.int32)).astype(np.uint32),
            "msp_widths": msp_sorted.get("widths", np.array([], dtype=np.int32)).astype(np.uint32),
        }

    # ------------------------------------------------------------------
    # HDF5 write helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _create_dataset(group: h5py.Group, name: str, data: np.ndarray) -> None:
        group.create_dataset(name, data=data, compression="gzip")

    def _write_chrom_data(
        self,
        sample_grp: h5py.Group,
        chrom: str,
        data: dict[str, np.ndarray | int],
    ) -> None:
        """Write one assembled chromosome to HDF5."""
        cg = sample_grp.create_group(chrom)
        cg.attrs["n_fibers"] = int(data["n_fibers"])

        self._create_dataset(cg, "fiber_id_table", data["fiber_id_table"])

        fmg = cg.create_group("fiber_metadata")
        self._create_dataset(fmg, "fiber_int_ids", data["fm_int_ids"])
        self._create_dataset(fmg, "starts", data["fib_starts"])
        self._create_dataset(fmg, "ends", data["fib_ends"])

        ng = cg.create_group("nucleosomes")
        self._create_dataset(ng, "fiber_int_ids", data["nuc_int_ids"])
        self._create_dataset(ng, "starts", data["nuc_starts"])
        self._create_dataset(ng, "ends", data["nuc_ends"])
        self._create_dataset(ng, "widths", data["nuc_widths"])
        self._create_dataset(ng, "linkers", data["nuc_linkers"])

        if len(data["mc_pos"]) > 0:
            mg = cg.create_group("5mC")
            self._create_dataset(mg, "fiber_int_ids", data["mc_int_ids"])
            self._create_dataset(mg, "positions", data["mc_pos"])
            self._create_dataset(mg, "probabilities", data["mc_prob"])
            self._create_dataset(mg, "is_methylated", data["mc_is_meth"])

        if len(data["hmc_pos"]) > 0:
            hg = cg.create_group("5hmC")
            self._create_dataset(hg, "fiber_int_ids", data["hmc_int_ids"])
            self._create_dataset(hg, "positions", data["hmc_pos"])
            self._create_dataset(hg, "probabilities", data["hmc_prob"])
            self._create_dataset(hg, "is_methylated", data["hmc_is_meth"])

        if len(data["ma_pos"]) > 0:
            ag = cg.create_group("6mA")
            self._create_dataset(ag, "fiber_int_ids", data["ma_int_ids"])
            self._create_dataset(ag, "positions", data["ma_pos"])

        if len(data["msp_starts"]) > 0:
            mg = cg.create_group("msp")
            self._create_dataset(mg, "fiber_int_ids", data["msp_int_ids"])
            self._create_dataset(mg, "starts", data["msp_starts"])
            self._create_dataset(mg, "ends", data["msp_ends"])
            self._create_dataset(mg, "widths", data["msp_widths"])

        sample_grp.file.flush()

    def _build_indices(
        self,
        sample_grp: h5py.Group,
        chrom: str,
        data: dict[str, np.ndarray | int],
    ) -> None:
        """Write per-layer fiber_int_id -> array slice indices."""
        ig = sample_grp[chrom].create_group("_indices")
        for layer, fid_arr in [
            ("nucleosomes", data["nuc_int_ids"]),
            ("5mC", data["mc_int_ids"]),
            ("5hmC", data["hmc_int_ids"]),
            ("6mA", data["ma_int_ids"]),
            ("msp", data["msp_int_ids"]),
        ]:
            if len(fid_arr) == 0:
                continue
            boundaries = np.concatenate([
                [0],
                np.where(fid_arr[1:] != fid_arr[:-1])[0] + 1,
                [len(fid_arr)],
            ])
            unique_fids = fid_arr[boundaries[:-1]]
            starts_idx = boundaries[:-1].astype(np.uint32)
            ends_idx = boundaries[1:].astype(np.uint32)
            slice_data = np.array(
                list(zip(unique_fids, starts_idx, ends_idx, strict=False)),
                dtype=schema.SLICE_INDEX_DTYPE,
            )
            self._create_dataset(ig, f"{layer}_slices", slice_data)

    # ------------------------------------------------------------------
    # Batch 2 stops here. Full sample orchestration and annotation loading
    # land in Batch 3.
    # ------------------------------------------------------------------

    def build(self) -> None:
        """Run the full build pipeline. To be implemented in Batch 3."""
        raise NotImplementedError(
            "FiberDatabaseBuilder.build() will be implemented in Batch 3."
        )
