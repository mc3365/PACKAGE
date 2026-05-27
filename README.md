# PACKAGE

[![Tests](https://github.com/YOUR_USERNAME/PACKAGE/actions/workflows/tests.yml/badge.svg)](https://github.com/YOUR_USERNAME/PACKAGE/actions/workflows/tests.yml)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**Single-molecule epigenomic data integration for long-read Fiber-seq.**

`PACKAGE` provides infrastructure for joint single-molecule analysis of DNA methylation,
nucleosome positioning, and chromatin accessibility on individual long reads from Oxford
Nanopore (ONT) and PacBio Fiber-seq experiments. It uses an HDF5-backed storage layer with
spatial indexing for genome-scale per-molecule queries.

## Status

🚧 **Pre-release (v0.1.0).** APIs are unstable; structure may change before the first
tagged release.

## What it does

`PACKAGE` integrates the outputs of established modification callers (`modkit`,
`fibertools-rs`, `pb-CpG-tools`) into a unified per-fiber data structure that supports
joint queries across modification types and platforms. Existing tools handle modification
extraction excellently but produce separate files in separate formats with no efficient
way to ask single-molecule questions across them. `PACKAGE` provides the missing
integration layer.

| Layer | Module | Description |
|-------|--------|-------------|
| **Extraction** | `PACKAGE.extract` | Wraps `modkit`, `fibertools-rs`, `pb-CpG-tools` with platform-aware defaults |
| **Storage** | `PACKAGE.db` | HDF5 schema with spatial indexing and memory-mapped access |
| **Visualization** | `PACKAGE.viz` | Per-molecule heatmap and related plots |
| **Benchmark** | `PACKAGE.benchmark` | Reproducible scripts for performance and cross-platform figures |

## Quick install

```bash
git clone https://github.com/YOUR_USERNAME/PACKAGE.git
cd PACKAGE
pip install -e ".[dev]"
```

## Quickstart

```bash
# Step 1: extract modifications from BAMs
PACKAGE extract --config configs/example.yaml --platform ont

# Step 2: build the HDF5 database
PACKAGE build --config configs/example.yaml

# Step 3: query
PACKAGE query --db fiber_database.h5 --region chr1:1000000-1100000
```

Or as a Python library:

```python
from PACKAGE.db import FiberDatabase
from PACKAGE.viz import single_molecule_heatmap

with FiberDatabase("fiber_database.h5") as db:
    fig = single_molecule_heatmap(db, "chr1", 1_000_000, 1_100_000, sample="d0")
```

See [`examples/`](examples/) for tutorial notebooks.

## Supported platforms

- **Oxford Nanopore (ONT)**: input from `modkit extract` + `ft extract` (fibertools-rs)
- **PacBio**: input from `pb-CpG-tools` (or fibertools-rs 5mC) + `ft extract`

Platform is selected via `--platform {ont,pacbio}` at extraction time. Downstream
storage, query, and analysis are platform-agnostic.

## Design philosophy

`PACKAGE` wraps rather than reimplements established modification callers. The novel
contribution is the integration layer: joint single-molecule storage, spatial indexing
across modification types, and the cross-platform abstraction. Analyses are demonstrated
through example notebooks rather than packaged as a library, leaving users free to
implement domain-specific analyses against the queryable data structure.

## Citation

If you use `PACKAGE` in published work, please cite:

> Cui, E. et al. (2026). *PACKAGE: a unified framework for single-molecule integration
> of long-read Fiber-seq data.* [Manuscript in preparation].

A `CITATION.cff` file is included for automated citation generation.

## License

MIT. See [LICENSE](LICENSE).

## Contact

Erika Cui — meiying.cui@yale.edu
