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

## HDF5 Layout

PACKAGE stores all samples and shared annotations in one HDF5 file. The file is
organized as three nested layers: file root, sample/chromosome groups, and flat
per-chromosome arrays. Within each chromosome, molecular features are sorted by
integer fiber ID. The `_indices` group maps each fiber ID to the row range for that
fiber in each feature array, so per-fiber lookups do not require scanning the full
chromosome.

![PACKAGE HDF5 data structure](figures/architecture/package_hdf5_database_structure.png)

## See also

- [Installation](installation.md) — getting started
- [Tutorial](tutorial.md) — building your first database
- [Parameter Reference](parameters.md) — command and YAML options
- [API Reference](api.md) — module documentation
