# PACKAGE

[![Tests](https://github.com/mc3365/PACKAGE/actions/workflows/tests.yml/badge.svg)](https://github.com/mc3365/PACKAGE/actions/workflows/tests.yml)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**Single-molecule epigenomic data integration for long-read Fiber-seq.**

`PACKAGE` provides infrastructure for joint single-molecule analysis of DNA methylation,
nucleosome positioning, and chromatin accessibility on individual long reads from Oxford
Nanopore (ONT) and PacBio Fiber-seq experiments. It uses an HDF5-backed storage layer
with spatial indexing for genome-scale per-molecule queries.

## Status

**Alpha (v0.4.0).** ONT extraction, HDF5 building, indexed queries, benchmarks, and
visualization have been validated end to end. PacBio extraction and normalization use
the same HDF5 schema for 5mC, 6mA, MSPs, nucleosomes, and optional FIRE accessibility
calls. The PacBio FIRE workflow can regenerate `Cov.bed`, build the enhancer-by-fiber
object, rank co-accessible pairs and stitched regions, and compare those outputs with
the legacy workflow. APIs may still change.

## What it does

`PACKAGE` integrates the outputs of established modification callers (`modkit`,
`fibertools-rs`) into a unified per-fiber data structure that supports joint queries
across modification types. Existing tools handle modification
extraction excellently but produce separate files in separate formats with no efficient
way to ask single-molecule questions across them. `PACKAGE` provides the missing
integration layer.

| Layer | Module | Description |
|-------|--------|-------------|
| **Extraction** | `PACKAGE.extract` | Wraps `modkit` and `fibertools-rs` with validated ONT defaults and PacBio normalization |
| **Storage** | `PACKAGE.db` | HDF5 schema with spatial indexing and memory-mapped access |
| **Visualization** | `PACKAGE.viz` | Per-molecule heatmap and related plots |
| **Co-accessibility** | `PACKAGE.analysis.coaccessibility` | PacBio FIRE coverage, object construction, and pair/region ranking |
| **Benchmark** | `PACKAGE.benchmark` | Query/storage benchmarks and legacy co-accessibility validation |

## Data Structure

PACKAGE writes one HDF5 file containing all samples and shared annotations. Each
sample is partitioned by chromosome, and each chromosome stores molecular features
as parallel flat arrays. Rows are sorted by integer fiber ID, and `_indices` maps
each fiber ID to its row range in every feature array. This keeps per-fiber lookups
fast without repeatedly scanning whole chromosomes.

![PACKAGE HDF5 data structure](docs/figures/architecture/package_hdf5_database_structure.png)

## Quick install

```bash
git clone https://github.com/mc3365/PACKAGE.git
cd PACKAGE
conda env create -f environments/ont.yml
conda activate package-ont
pip install -e ".[dev,viz]"
```

## Quickstart

```bash
# Step 1: extract modifications from BAMs
PACKAGE extract --platform ont --config configs/my_ont.yaml --samples d0

# Step 2: build the HDF5 database and spatial-index sidecar
PACKAGE build --config configs/my_ont.yaml --build-index

# Step 3: query one region
PACKAGE query \
  --db /path/to/output/fiber_database.h5 \
  --region chr1:1000000-1100000 \
  --sample d0
```

Or as a Python library:

```python
from PACKAGE.db import FiberDatabase
from PACKAGE.viz import single_molecule_heatmap

with FiberDatabase("fiber_database.h5") as db:
    fig = single_molecule_heatmap(db, "chr1", 1_000_000, 1_100_000, sample="d0")
    fig.savefig("ont_region.png", dpi=200)
```

For PacBio data, use `PACKAGE extract --platform pacbio` with a PacBio config, or
normalize existing `ft extract` outputs before building the database. See the
[tutorial](docs/tutorial.md#pacbio-workflow) for the tested end-to-end workflow.

PacBio co-accessibility support starts from a FIRE-annotated BAM or the
`acc.model.results.bed` generated from it by `ft fire --extract`. The full upstream
FIRE snakemake/modeling workflow is not bundled; for exact reproduction of legacy
co-accessibility outputs, use the same extracted FIRE BED and fibertools version
provenance as the original run.

Start with the [tutorial](docs/tutorial.md), then use the
[parameter reference](docs/parameters.md) for every YAML and CLI option. The
[benchmark and validation guide](docs/benchmark.md) explains outputs and figures,
and the [step-by-step ONT workflow](STEP_BY_STEP_GITHUB.md) provides a focused ONT
walkthrough.

## Supported platforms

- **Oxford Nanopore (ONT)**: validated with `modkit extract` and `ft extract`
- **PacBio**: validated for fibertools-derived 5mC, 6mA, MSP, nucleosome, and
  optional `ft fire --extract` accessibility-call layers after normalization into the
  shared HDF5 schema; FIRE co-accessibility coverage, object construction, ranking,
  and legacy-output validation are also available

The package stores the FIRE/accessibility input layer needed by the legacy
co-accessibility workflow and can regenerate the legacy `Cov.bed`,
enhancer-by-fiber JSON object, and ranked constituent-element pair tables from that
layer. PACKAGE does not run the upstream FIRE model that creates a FIRE-annotated
BAM. More specialized downstream biological interpretation remains workflow-specific.

## Design philosophy

`PACKAGE` wraps rather than reimplements established modification callers. The novel
contribution is the integration layer: joint single-molecule storage, spatial indexing
across modification types, and the cross-platform abstraction. Reusable operations
such as annotation summaries, single-molecule visualization, and FIRE
co-accessibility ranking are packaged, while study-specific biological interpretation
remains downstream.

## Citation

If you use `PACKAGE` in published work, please cite:

> Cui, M. et al. (2026). *PACKAGE: a unified framework for single-molecule integration
> of long-read Fiber-seq data.* [Manuscript in preparation].

A `CITATION.cff` file is included for automated citation generation.

## License

MIT. See [LICENSE](LICENSE).

## Contact

Meiying Cui — meiying.cui@yale.edu
