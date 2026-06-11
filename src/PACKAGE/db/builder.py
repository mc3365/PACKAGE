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

from tqdm import tqdm

from PACKAGE.config import Config, SampleConfig
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
    # Batch 1 stops here. The class-level skeleton continues below as stubs
    # so callers see clear NotImplementedError messages rather than
    # AttributeError.
    # ------------------------------------------------------------------

    def build(self) -> None:
        """Run the full build pipeline. To be implemented in Batch 3."""
        raise NotImplementedError(
            "FiberDatabaseBuilder.build() will be implemented in Batch 3."
        )
