# API Reference

API docs are auto-generated from docstrings via `mkdocstrings`.

The database API is built around a chromosome-partitioned HDF5 layout. Samples and
shared annotations live in one file, while each sample/chromosome group stores flat
feature arrays plus `_indices` tables that map integer fiber IDs to row ranges.
The tutorial includes a [query walkthrough](tutorial.md#how-a-query-uses-the-database)
showing how this layout is used for annotation-centered queries.

## Database

::: PACKAGE.db.database
    options:
      show_source: false
      members:
        - FiberDatabase

## Configuration

::: PACKAGE.config
    options:
      show_source: false

## Extraction

::: PACKAGE.extract.ont
    options:
      show_source: false

::: PACKAGE.extract.pacbio
    options:
      show_source: false

## Co-Accessibility

::: PACKAGE.analysis.coaccessibility
    options:
      show_source: false

## Benchmarks

::: PACKAGE.benchmark.annotation
    options:
      show_source: false

::: PACKAGE.benchmark.coaccessibility
    options:
      show_source: false

## Visualization

::: PACKAGE.viz.single_molecule
    options:
      show_source: false

## Schema

::: PACKAGE.db.schema
    options:
      show_source: true
