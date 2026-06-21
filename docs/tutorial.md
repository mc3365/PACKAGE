# ONT Tutorial

## Prepare the configuration

Copy `configs/ont_template.yaml`, replace every `/path/to` value, and keep the 5mC
and 5hmC layer paths identical. The shared modkit table contains both call types.
The configured methylation threshold is `0.5`; binarization happens during the HDF5
build while raw probabilities are retained.

The BAM must be coordinate sorted, indexed, aligned to the configured reference, and
contain MM/ML tags. PACKAGE checks these requirements before extraction.

## Extract and build

```bash
PACKAGE extract --platform ont --config configs/my_ont.yaml
PACKAGE build --config configs/my_ont.yaml
PACKAGE info /path/to/output/fiber_database.h5
```

Existing extraction outputs are reused when `overwrite: false`. Each extraction writes
`PACKAGE_manifest_<sample>.json` with input QC, tool versions, commands, and outputs.

## Query a region

```bash
PACKAGE query \
  --db /path/to/output/fiber_database.h5 \
  --region chr1:3000000-5000000 \
  --sample sample1
```

## Plot individual molecules

```bash
python examples/ont_region_plot.py \
  --db /path/to/output/fiber_database.h5 \
  --region chr1:3000000-3050000 \
  --sample sample1 \
  --layer 5mC \
  --out ont_region.png
```

The figure shows each fiber as a row, with nucleosomes, MSP intervals, and base-level
modification calls. Limit plots to focused regions and use `--max-fibers` for legibility.
