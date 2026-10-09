# API Reference

API docs are auto-generated from docstrings via `mkdocstrings`.

The database API is built around a chromosome-partitioned HDF5 layout. Samples and
shared annotations live in one file, while each sample/chromosome group stores flat
feature arrays plus `_indices` tables that map integer fiber IDs to row ranges.
The tutorial includes a [query walkthrough](tutorial.md#how-a-query-uses-the-database)
showing how this layout is used for annotation-centered queries.

## Database

::: mei_fiber.db.database
    options:
      show_source: false
      members:
        - FiberDatabase

## Configuration

::: mei_fiber.config
    options:
      show_source: false

## Extraction

::: mei_fiber.extract.ont
    options:
      show_source: false

::: mei_fiber.extract.pacbio
    options:
      show_source: false

## Co-Accessibility

::: mei_fiber.analysis.coaccessibility
    options:
      show_source: false

## Benchmarks

::: mei_fiber.benchmark.annotation
    options:
      show_source: false

::: mei_fiber.benchmark.coaccessibility
    options:
      show_source: false

## Visualization

::: mei_fiber.viz.single_molecule
    options:
      show_source: false

## Schema

::: mei_fiber.db.schema
    options:
      show_source: true
