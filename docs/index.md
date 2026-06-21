# PACKAGE

Single-molecule epigenomic analysis for ONT long-read Fiber-seq data.

## Overview

`PACKAGE` provides a unified analytical framework for joint analysis of:

- DNA methylation (5mC, 5hmC)
- 6mA-mediated chromatin accessibility
- Nucleosome positioning

...on individual Oxford Nanopore (ONT) long reads.

The package uses an HDF5-backed storage layer with spatial indexing for genome-scale
per-molecule queries.

ONT extraction, database building, regional queries, and visualization are available.
PacBio extraction is planned for a later release.

## See also

- [Installation](installation.md) — getting started
- [Tutorial](tutorial.md) — building your first database
- [API Reference](api.md) — module documentation
