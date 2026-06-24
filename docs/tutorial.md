# ONT Tutorial

## Prepare the configuration

Copy `configs/ont_template.yaml`, replace every `/path/to` value, and keep the 5mC
and 5hmC layer paths identical. The shared modkit table contains both call types.
The configured methylation threshold is `0.5`; binarization happens during the HDF5
build while raw probabilities are retained.

The BAM must be coordinate sorted, indexed, aligned to the configured reference, and
contain MM/ML tags. PACKAGE checks these requirements before extraction.

For all YAML fields and command-line options, see the
[Parameter Reference](parameters.md).

## Extract and build

```bash
PACKAGE extract --platform ont --config configs/my_ont.yaml
PACKAGE build --config configs/my_ont.yaml
PACKAGE info /path/to/output/fiber_database.h5
```

Existing extraction outputs are reused when `overwrite: false`. Each extraction writes
`PACKAGE_manifest_<sample>.json` with input QC, tool versions, commands, and outputs.
With `build.build_spatial_index: true`, the build also writes
`fiber_database.index.pkl` beside the HDF5 file. Keep the sidecar with the database;
it is used for fast coordinate queries.

## Query a region

```bash
PACKAGE query \
  --db /path/to/output/fiber_database.h5 \
  --region chr1:3000000-5000000 \
  --sample sample1
```

## Visualization Examples

Install the optional plotting dependencies before running these examples:

```bash
pip install -e ".[viz]"
```

The examples below use an existing PACKAGE HDF5 database. They are intended as
small, inspectable outputs rather than final manuscript layouts.

### ECDF Summary Plots

Global ECDF plots summarize per-fiber feature fractions for one or more samples:

```bash
python examples/ont_feature_ecdf.py \
  --db /path/to/output/fiber_database.h5 \
  --samples d0 d4 \
  --outdir figures/ecdf
```

This writes:

- `global_feature_fractions.csv`
- `ecdf_coverage_fraction.png/.pdf`
- `ecdf_base_specific_fraction.png/.pdf`

Each row of `global_feature_fractions.csv` is one fiber. Nucleosome and MSP
coverage are measured as covered base pairs divided by fiber length. 5mC, 5hmC,
and 6mA coverage fractions are call counts divided by fiber length. The
base-specific panel reports 5mC and 5hmC fractions among CpG calls using the
stored binary calls; the 6mA panel is a density-style proxy because the database
stores called 6mA positions, not every adenine.

The ECDF y-axis is the fraction of fibers at or below each x-axis value, so left
or right shifts between samples indicate global differences in per-fiber feature
burden.

![Per-fiber feature coverage ECDF](figures/visualization/ecdf/ecdf_coverage_fraction.png)

![Per-fiber base-call fraction ECDF](figures/visualization/ecdf/ecdf_base_specific_fraction.png)

Use `--max-fibers-per-chrom` for a fast smoke test:

```bash
python examples/ont_feature_ecdf.py \
  --db /path/to/output/fiber_database.h5 \
  --samples d0 d4 \
  --max-fibers-per-chrom 5000 \
  --outdir figures/ecdf_smoke
```

### Centered Heatmap and Metaplot

For fixed-window center/TSS plots, use the centered workflow. Pass the 9-column
annotation BED when promoter or gene-body strand orientation should be respected:

```bash
python examples/ont_centered_heatmap.py \
  --db /path/to/output/fiber_database.h5 \
  --bed /path/to/master_annotations_v4.uniqueID.bed \
  --annotation CGI \
  --samples d0 d4 \
  --methylation-display binned \
  --outdir figures/CGI_center
```

This writes one heatmap per sample and one multi-sample metaplot:

- `<annotation>_<sample>_center_heatmap.png/.pdf`
- `<annotation>_center_metaplot.png/.pdf`
- `<annotation>_<sample>_nuc.npz`
- `<annotation>_<sample>_met.npz`
- `<annotation>_<sample>_met_binned.npz`
- `<annotation>_<sample>_metaplot.csv`
- `region_summary.csv`

For CGI and other non-directional regions, the window is centered at the
annotation midpoint. For promoter-like annotations, the window is centered at
the strand-aware TSS when the BED contains strand information. Minus-strand
directional regions are flipped before averaging so upstream remains on the left.

The heatmap rows are retained annotation regions after averaging across fibers
that span each region. The left heatmap shows nucleosome occupancy, the right
heatmap shows 5mC, and the side bars show row means. With
`--methylation-display binned`, 5mC is shown in 50 bp bins, which gives a smoother
matrix for sparse CpG calls.

![CGI centered heatmap](figures/visualization/heatmap/CGI_d0_center_heatmap.png)

The metaplot averages the retained region-level rows into one profile per sample
and signal. It is useful for checking whether the expected centered structure is
recovered across many annotations.

![CGI centered metaplot](figures/visualization/heatmap/CGI_center_metaplot.png)

For a faster test, cap the number of regions:

```bash
python examples/ont_centered_heatmap.py \
  --db /path/to/output/fiber_database.h5 \
  --bed /path/to/master_annotations_v4.uniqueID.bed \
  --annotation CGI \
  --samples d0 d4 \
  --methylation-display binned \
  --max-regions 50 \
  --outdir figures/CGI_center_smoke
```

For body-normalized annotation plots rather than fixed center/TSS windows, use
`examples/ont_annotation_heatmap.py`.

### Single-Molecule Region Plot

```bash
python examples/ont_region_plot.py \
  --db /path/to/output/fiber_database.h5 \
  --region chr1:3000000-3050000 \
  --sample sample1 \
  --layer 5mC \
  --out ont_region.png
```

The figure shows each fiber as one row. The thin black line marks the displayed
span of the fiber, gray blocks mark nucleosomes, teal outlines mark MSPs, and
vertical ticks mark the selected modification layer. For 5mC and 5hmC, red ticks
are modified calls and pale gray ticks are unmodified calls.

```bash
python examples/ont_region_plot.py \
  --db /path/to/output/fiber_database.h5 \
  --region chr1:3000000-3050000 \
  --sample d0 \
  --layer 5mC \
  --tracks modification \
  --hide-unmodified \
  --max-fibers 100 \
  --out figures/single_molecule/d0_chr1_5mC_clean.png
```

Use track modes to make the plot less crowded:

```bash
python examples/ont_region_plot.py \
  --db /path/to/output/fiber_database.h5 \
  --region chr1:3000000-3050000 \
  --sample d0 \
  --tracks chromatin \
  --max-fibers 100 \
  --out figures/single_molecule/d0_chr1_chromatin_only.png
```

Available track modes are:

- `full`: nucleosomes, MSPs, and the selected modification layer;
- `modification`: selected modification layer only;
- `chromatin`: nucleosomes and MSPs only.

Limit single-molecule plots to focused windows and use `--max-fibers` for
legibility.

![Single-molecule 5mC view](figures/visualization/single_molecule/d0_chr1_3000000_3050000_5mC_clean.png)

![Single-molecule chromatin-only view](figures/visualization/single_molecule/d0_chr1_chromatin_only.png)
