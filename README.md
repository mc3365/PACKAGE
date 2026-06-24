# PACKAGE

[![Tests](https://github.com/mc3365/PACKAGE/actions/workflows/tests.yml/badge.svg)](https://github.com/mc3365/PACKAGE/actions/workflows/tests.yml)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**Single-molecule epigenomic data integration for long-read Fiber-seq.**

`PACKAGE` provides infrastructure for joint single-molecule analysis of DNA methylation,
nucleosome positioning, and chromatin accessibility on individual long reads from Oxford
Nanopore (ONT) experiments. It uses an HDF5-backed storage layer with
spatial indexing for genome-scale per-molecule queries.

## Status

**ONT alpha (v0.4.0).** ONT extraction, HDF5 building, regional queries, and
single-molecule visualization have been validated end to end. APIs may still change.

## What it does

`PACKAGE` integrates the outputs of established modification callers (`modkit`,
`fibertools-rs`) into a unified per-fiber data structure that supports joint queries
across modification types. Existing tools handle modification
extraction excellently but produce separate files in separate formats with no efficient
way to ask single-molecule questions across them. `PACKAGE` provides the missing
integration layer.

| Layer | Module | Description |
|-------|--------|-------------|
| **Extraction** | `PACKAGE.extract` | Wraps `modkit` and `fibertools-rs` with validated ONT defaults |
| **Storage** | `PACKAGE.db` | HDF5 schema with spatial indexing and memory-mapped access |
| **Visualization** | `PACKAGE.viz` | Per-molecule heatmap and related plots |
| **Benchmark** | `PACKAGE.benchmark` | Reproducible scripts for performance and cross-platform figures |

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

See the [step-by-step ONT workflow](STEP_BY_STEP_GITHUB.md),
the [parameter reference](docs/parameters.md), and [`examples/`](examples/)
for complete usage examples.

## Supported platforms

- **Oxford Nanopore (ONT)**: validated with `modkit extract` and `ft extract`
- **PacBio**: planned; its extraction adapter is not implemented in this release

## Design philosophy

`PACKAGE` wraps rather than reimplements established modification callers. The novel
contribution is the integration layer: joint single-molecule storage, spatial indexing
across modification types, and the cross-platform abstraction. Analyses are demonstrated
through example notebooks rather than packaged as a library, leaving users free to
implement domain-specific analyses against the queryable data structure.

## Citation

If you use `PACKAGE` in published work, please cite:

> Cui, M. et al. (2026). *PACKAGE: a unified framework for single-molecule integration
> of long-read Fiber-seq data.* [Manuscript in preparation].

A `CITATION.cff` file is included for automated citation generation.

## License

MIT. See [LICENSE](LICENSE).

## Contact

Meiying Cui — meiying.cui@yale.edu
