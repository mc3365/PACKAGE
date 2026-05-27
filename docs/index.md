# PACKAGE

Single-molecule epigenomic analysis for long-read Fiber-seq data.

## Overview

`PACKAGE` provides a unified analytical framework for joint analysis of:

- DNA methylation (5mC, 5hmC)
- 6mA-mediated chromatin accessibility
- Nucleosome positioning

...on individual long reads from Oxford Nanopore (ONT) and PacBio Fiber-seq experiments.

The package uses an HDF5-backed storage layer with spatial indexing for genome-scale
per-molecule queries.

## See also

- [Installation](installation.md) — getting started
- [Tutorial](tutorial.md) — building your first database
- [API Reference](api.md) — module documentation
