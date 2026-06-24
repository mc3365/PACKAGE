# Parameter Reference

This page lists the user-facing parameters for the ONT workflow. It is organized by
workflow stage: extraction, HDF5 build, query, benchmark, and plotting.

PACKAGE currently validates the ONT workflow. PacBio extraction is planned for a
future release.

## YAML Configuration

Start from `configs/ont_template.yaml` and replace every `/path/to` value.

| Field | Required | Meaning |
| --- | --- | --- |
| `output_dir` | yes | Directory where the HDF5 database is written. |
| `output_file` | yes | HDF5 filename created by `PACKAGE build`. |
| `reference` | yes for extraction | Reference FASTA used by `modkit extract full --reference`. It should match the BAM genome assembly. |
| `samples[].name` | yes | Short sample name. This becomes the top-level sample name in the HDF5 database. |
| `samples[].bam` | yes for extraction | Coordinate-sorted, indexed, aligned BAM with MM/ML modification tags. |
| `samples[].layers.nucleosomes` | yes for build | Flattened nucleosome CSV produced by `PACKAGE extract`. |
| `samples[].layers.5mC` | optional | Modkit `extract full` TSV or TSV.GZ. Usually the same file as `5hmC`. |
| `samples[].layers.5hmC` | optional | Modkit `extract full` TSV or TSV.GZ. Usually the same file as `5mC`. |
| `samples[].layers.6mA` | optional | Fibertools `ft extract --m6a` BED output. |
| `samples[].layers.msp` | optional | Fibertools `ft extract --msp` BED output. |
| `annotations` | optional | Named BED4 annotation files loaded into the HDF5 database. |
| `parameters.methylation_threshold` | yes | Probability cutoff used to create binary 5mC/5hmC calls during HDF5 build. The validated ONT value is `0.5`. Raw probabilities are also retained. |
| `parameters.min_msp_size` | yes | Minimum MSP interval size kept during HDF5 build. |
| `build.build_spatial_index` | recommended | Whether `PACKAGE build` should also create `<database>.index.pkl` for fast coordinate queries. |

### Extraction YAML

| Field | Default in template | Meaning |
| --- | --- | --- |
| `extraction.threads` | `16` | Threads passed to `modkit`. |
| `extraction.overwrite` | `false` | When `false`, existing extraction outputs are reused. Set `true` to regenerate them. |
| `extraction.keep_nucleosome_bed` | `true` | Keep the intermediate fibertools BED12 nucleosome file beside the flattened CSV. |
| `extraction.validate_bam` | `true` | Check BAM readability, index, and sampled MM/ML tags before extraction. |
| `extraction.write_manifest` | `true` | Write `PACKAGE_manifest_<sample>.json` with tool versions, commands, QC, and output sizes. |
| `extraction.modkit_executable` | `modkit` | Executable name or absolute path for modkit. |
| `extraction.ft_executable` | `ft` | Executable name or absolute path for fibertools-rs. |
| `extraction.samtools_executable` | `samtools` | Executable name or absolute path for samtools. |

Use absolute executable paths on HPC systems when the Python environment and ONT
tool environment are different.

## Extraction

```bash
PACKAGE extract \
  --platform ont \
  --config configs/my_ont.yaml \
  --samples d0
```

| Option | Required | Meaning |
| --- | --- | --- |
| `--platform ont` | yes | Selects the ONT extractor. `pacbio` is reserved for future support. |
| `--config` | yes | YAML configuration file. |
| `--samples` | no | Sample name to extract. Repeat the flag for multiple samples. If omitted, all samples in the YAML are processed. |

Outputs per sample:

- `<sample>_raw_mods.tsv.gz`: modkit table containing 5mC and 5hmC rows.
- `<sample>_6ma.bed`: fibertools 6mA intervals.
- `<sample>_msp.bed`: fibertools MSP intervals.
- `<sample>_nuc_features.bed12`: raw fibertools nucleosome BED12 when kept.
- `<sample>_nuc_features.csv`: flattened nucleosome table used by the builder.
- `PACKAGE_manifest_<sample>.json`: extraction QC and provenance when enabled.

## Build

```bash
PACKAGE build \
  --config configs/my_ont.yaml \
  --samples d0 \
  --build-index
```

| Option | Required | Meaning |
| --- | --- | --- |
| `--config` | yes | YAML configuration file describing extracted layers and output database path. |
| `--samples` | no | Sample name to build. Repeat for multiple samples. If omitted, all samples in the YAML are built. |
| `--build-index` | no | Build the cached spatial-index sidecar after the HDF5 file is created. |
| `--no-build-index` | no | Skip spatial-index creation even if enabled in the YAML. |

The build writes `output_dir/output_file`. If spatial indexing is enabled, it also
writes `<output_file>.index.pkl` beside the HDF5 file. Keep the sidecar with the HDF5
database for fast region queries.

## Info And Region Query

```bash
PACKAGE info /path/to/fiber_database.h5
```

`PACKAGE info` reports version, genome, samples, annotation count, spatial-index
status, and per-sample fiber counts.

```bash
PACKAGE query \
  --db /path/to/fiber_database.h5 \
  --region chr1:3000000-5000000 \
  --sample d0 \
  --out d0_chr1_fibers.txt
```

| Option | Required | Meaning |
| --- | --- | --- |
| `--db` | yes | HDF5 database path. |
| `--region` | yes | Genomic interval in `chrom:start-end` form. |
| `--sample` | no | Restrict query to one sample. If omitted, all samples are searched. |
| `--out` | no | Write overlapping fiber IDs to a text file. If omitted, only the count is printed. |

## Random-Region Benchmark

```bash
python -m PACKAGE.benchmark.ont \
  --db /path/to/fiber_database.h5 \
  --config configs/my_ont.yaml \
  --sample d0 \
  --outdir benchmark/d0 \
  --region-sizes-kb 1 10 100 1000 \
  --n-regions 25 \
  --repeats 3 \
  --seed 3365 \
  --query-modes spatial_index \
  --build-index
```

| Option | Default | Meaning |
| --- | --- | --- |
| `--db` | required | HDF5 database to benchmark. |
| `--outdir` | required | Directory for CSV, JSON, PNG, and PDF outputs. |
| `--sample` | first sample | Sample used to choose and query regions. |
| `--config` | none | YAML file used to report intermediate extraction file sizes. |
| `--region-sizes-kb` | `1 10 100 1000` | Window sizes tested around sampled fibers. |
| `--n-regions` | `25` | Number of regions per window size. |
| `--repeats` | `3` | Timed repeats per selected region. |
| `--seed` | `3365` | Random seed for reproducible region selection. |
| `--query-modes` | `array_scan spatial_index` | Query implementations to time. Use `spatial_index` alone for full production databases. |
| `--build-index` | off | Build the spatial-index sidecar before timing. |

The random benchmark can be slow if `array_scan` is included on a full database. For
routine production checks, use `--query-modes spatial_index`.

## Annotation Benchmark

```bash
python -m PACKAGE.benchmark.annotation \
  --db /path/to/fiber_database.h5 \
  --sample d0 \
  --annotation CGI \
  --layers nucleosomes 5mC 5hmC msp \
  --max-regions 50 \
  --repeats 3 \
  --warmups 1 \
  --outdir benchmark/d0_CGI
```

| Option | Default | Meaning |
| --- | --- | --- |
| `--db` | required | HDF5 database to benchmark. |
| `--outdir` | required | Directory for annotation benchmark outputs. |
| `--sample` | required | Sample queried. |
| `--annotation` | required | Annotation class to query. Repeat for multiple classes. |
| `--layers` | `nucleosomes 5mC` | Layers passed to `FiberDatabase.query_annotation_fast`. |
| `--max-regions` | all regions | Limit the number of annotation regions. Useful for smoke tests. |
| `--repeats` | `3` | Timed repeats. |
| `--warmups` | `1` | Untimed warm-up queries before timing. |

This benchmark is usually easiest to interpret for method figures because it times a
real biological query pattern.

## Plotting

### Single-Molecule Region Plot

```bash
python examples/ont_region_plot.py \
  --db /path/to/fiber_database.h5 \
  --region chr1:3000000-3050000 \
  --sample d0 \
  --layer 5mC \
  --tracks modification \
  --max-fibers 100 \
  --hide-unmodified \
  --out figures/d0_chr1_region.png
```

| Option | Default | Meaning |
| --- | --- | --- |
| `--db` | required | HDF5 database. |
| `--region` | required | Focused interval in `chrom:start-end` form. |
| `--sample` | first sample | Sample to plot. |
| `--layer` | `5mC` | Modification layer: `5mC`, `5hmC`, or `6mA`. |
| `--tracks` | `full` | Track set: `full`, `modification`, or `chromatin`. |
| `--max-fibers` | `100` | Maximum number of fiber rows shown. |
| `--hide-unmodified` | off | Hide unmodified 5mC/5hmC calls for a cleaner modification-focused view. |
| `--out` | `ont_region.png` | Output figure path. |

Use `--tracks modification` for 5mC/5hmC/6mA-only plots, `--tracks chromatin` for
nucleosome/MSP-only plots, and `--tracks full` for all tracks together.

### Global ECDF Plots

```bash
python examples/ont_feature_ecdf.py \
  --db /path/to/fiber_database.h5 \
  --samples d0 d4 \
  --max-fibers-per-chrom 50000 \
  --outdir figures/ecdf
```

| Option | Default | Meaning |
| --- | --- | --- |
| `--db` | required | HDF5 database. |
| `--outdir` | required | Output directory. |
| `--samples` | all samples | One or more samples to include. |
| `--max-fibers-per-chrom` | all fibers | Optional downsampling cap per chromosome for faster exploratory plots. |

### Normalized Annotation Heatmap

```bash
python examples/ont_annotation_heatmap.py \
  --db /path/to/fiber_database.h5 \
  --annotation CGI \
  --sample d0 \
  --max-regions 50 \
  --min-fibers 10 \
  --outdir figures/CGI_d0
```

| Option | Default | Meaning |
| --- | --- | --- |
| `--db` | required | HDF5 database. |
| `--outdir` | required | Output directory for matrices and figures. |
| `--annotation` | required | Annotation class stored in the database. |
| `--sample` | required | Sample to plot. |
| `--max-regions` | all regions | Limit regions for smoke tests. Omit for full annotation plots. |
| `--min-fibers` | `10` | Minimum spanning fibers required for a region. |

This workflow rescales each annotation body to fixed bins and includes flanking bins
for the metaplot.

### Centered Annotation Heatmap And Metaplot

```bash
python examples/ont_centered_heatmap.py \
  --db /path/to/fiber_database.h5 \
  --bed /path/to/master_annotations_v4.uniqueID.bed \
  --annotation CGI \
  --samples d0 d4 \
  --methylation-display binned \
  --min-fibers 10 \
  --outdir figures/CGI_center
```

| Option | Default | Meaning |
| --- | --- | --- |
| `--db` | required | HDF5 database. |
| `--outdir` | required | Output directory for centered matrices and figures. |
| `--annotation` | required | Annotation class. |
| `--samples` | required | One or more samples. |
| `--bed` | none | Optional 9-column annotation BED used for strand-aware promoter/TSS orientation. |
| `--max-regions` | all regions | Limit regions for smoke tests. Omit for full annotation plots. |
| `--min-fibers` | `10` | Minimum spanning fibers required for a region. |
| `--methylation-display` | `binned` | Heatmap methylation matrix: `binned` for 50 bp bins or `smoothed` for base-resolution smoothed values. |

For CGI, the centered workflow uses a fixed +/-2.5 kb window around the annotation
midpoint. Promoter classes use a strand-aware TSS center when a 9-column BED is
provided.
