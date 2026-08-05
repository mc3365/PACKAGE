# Parameter Reference

This page lists user-facing parameters for ONT and PacBio Phase 1 workflows. It is
organized by workflow stage: extraction, HDF5 build, query, benchmark, and plotting.

ONT is validated end to end. PacBio Phase 1 normalizes fibertools outputs into the
same HDF5 schema for nucleosomes, 5mC, 6mA, and MSPs. PacBio 5hmC and FIRE
score/co-accessibility HDF5 layers are not yet built.

## YAML Configuration

Start from `configs/ont_template.yaml` or `configs/pacbio_template.yaml` and replace
every `/path/to` value.

| Field | Required | Meaning |
| --- | --- | --- |
| `output_dir` | yes | Directory where the HDF5 database is written. |
| `output_file` | yes | HDF5 filename created by `PACKAGE build`. |
| `reference` | yes for extraction | Reference FASTA used by ONT `modkit extract full --reference`. It should match the BAM genome assembly. |
| `samples[].name` | yes | Short sample name. This becomes the top-level sample name in the HDF5 database. |
| `samples[].bam` | yes for extraction | Coordinate-sorted, indexed, aligned BAM. ONT BAMs must have MM/ML modification tags; PacBio BAMs should be fibertools-compatible Fiber-seq or FIRE BAMs. |
| `samples[].layers.nucleosomes` | yes for build | Flattened nucleosome CSV produced by `PACKAGE extract`. |
| `samples[].layers.5mC` | optional | ONT: modkit `extract full` TSV or TSV.GZ, usually the same file as `5hmC`. PacBio: PACKAGE-normalized modkit-like TSV/TSV.GZ created from `ft extract --all`. |
| `samples[].layers.5hmC` | optional | ONT only: modkit `extract full` TSV or TSV.GZ. Usually the same file as `5mC`. |
| `samples[].layers.6mA` | optional | Fibertools `ft extract --m6a` BED output. PacBio extraction normalizes this to builder-compatible BED12. |
| `samples[].layers.msp` | optional | Fibertools `ft extract --msp` BED output. PacBio extraction normalizes this to builder-compatible BED12. |
| `annotations` | optional | Named BED annotation files loaded into the HDF5 database during build. Rebuild the database to add annotations if they were omitted. |
| `parameters.methylation_threshold` | yes | Probability cutoff used to create binary 5mC/5hmC calls during HDF5 build. The validated ONT value is `0.5`. Raw probabilities are also retained. |
| `parameters.min_msp_size` | yes | Minimum MSP interval size kept during HDF5 build. |
| `build.build_spatial_index` | recommended | Whether `PACKAGE build` should also create `<database>.index.pkl` for fast coordinate queries. |

### Extraction YAML

| Field | Default in template | Meaning |
| --- | --- | --- |
| `extraction.threads` | `16` | Threads used by platform extraction tools when supported. |
| `extraction.overwrite` | `false` | When `false`, existing extraction outputs are reused. Set `true` to regenerate them. |
| `extraction.keep_nucleosome_bed` | `true` | Keep the intermediate fibertools BED12 nucleosome file beside the flattened CSV. |
| `extraction.validate_bam` | `true` | Check BAM readability, index, and sampled MM/ML tags before extraction. |
| `extraction.write_manifest` | `true` | Write `PACKAGE_manifest_<sample>.json` with tool versions, commands, QC, and output sizes. |
| `extraction.modkit_executable` | `modkit` | Executable name or absolute path for modkit. |
| `extraction.ft_executable` | `ft` | Executable name or absolute path for fibertools-rs. |
| `extraction.samtools_executable` | `samtools` | Executable name or absolute path for samtools. |

Use absolute executable paths on HPC systems when the Python environment and
long-read tool environment are different.

## Extraction

```bash
PACKAGE extract \
  --platform ont \
  --config configs/my_ont.yaml \
  --samples d0
```

| Option | Required | Meaning |
| --- | --- | --- |
| `--platform ont/pacbio` | yes | Selects the platform-specific extractor. |
| `--config` | yes | YAML configuration file. |
| `--samples` | no | Sample name to extract. Repeat the flag for multiple samples. If omitted, all samples in the YAML are processed. |

Outputs per sample:

ONT outputs:

- `<sample>_raw_mods.tsv.gz`: modkit table containing 5mC and 5hmC rows.
- `<sample>_6ma.bed`: fibertools 6mA intervals.
- `<sample>_msp.bed`: fibertools MSP intervals.
- `<sample>_nuc_features.bed12`: raw fibertools nucleosome BED12 when kept.
- `<sample>_nuc_features.csv`: flattened nucleosome table used by the builder.
- `PACKAGE_manifest_<sample>.json`: extraction QC and provenance when enabled.

PacBio outputs:

- `<sample>_5mc_for_PACKAGE.tsv.gz`: modkit-like 5mC TSV created from `ft extract --all`.
- `<sample>_6ma_for_PACKAGE.bed`: normalized BED12 6mA file.
- `<sample>_msp_for_PACKAGE.bed`: normalized BED12 MSP file.
- `<sample>_nuc_features.csv`: flattened nucleosome table used by the builder.
- raw sidecar files such as `<sample>.fibertools_all.tsv.gz` and `<sample>.raw.bed` may be kept for inspection.
- `PACKAGE_manifest_<sample>.json`: extraction QC, conversion QC, and provenance when enabled.

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

Annotations and molecular layers are written during the build. PACKAGE does not yet
provide a command to append annotations or FIRE/co-accessibility layers to an existing
HDF5 file in place; rebuild from the same intermediate files when adding those inputs.

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

## Benchmark Summary Figures

```bash
python examples/ont_benchmark_summary.py \
  --storage-summary benchmark/d0/storage_summary.json \
  --query-results benchmark/d0/query_speed_results.csv \
  --annotation-summary benchmark/annotation_d0/annotation_query_summary.csv \
  --bam /path/to/d0.output.fiberseq.bam \
  --bam /path/to/d4.output.fiberseq.bam \
  --outdir figures/benchmark_summary
```

| Option | Default | Meaning |
| --- | --- | --- |
| `--storage-summary` | required | `storage_summary.json` from `PACKAGE.benchmark.ont`. |
| `--query-results` | none | Optional `query_speed_results.csv` for random-window scatter plots. |
| `--annotation-summary` | none | Optional `annotation_query_summary.csv` for annotation timing bars. |
| `--bam` | none | Optional aligned BAM path. Repeat for multiple BAMs. |
| `--outdir` | required | Output directory for summary figures. |
| `--outlier-fibers` | `10000` | Returned-fiber threshold used to mark random-query outliers. |

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
| `--sample-colors` | automatic | Optional sample colors as `SAMPLE=COLOR` entries, for example `WT=#3b6fb6 KO=#c74f46`. |

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
| `--bed` | none | Optional matching 9-column annotation BED used for strand-aware promoter/TSS orientation. Use the same annotation file used to build the database so chromosome/start/end/region ID records match. |
| `--max-regions` | all regions | Limit regions for smoke tests. Omit for full annotation plots. |
| `--min-fibers` | `10` | Minimum spanning fibers required for a region. |
| `--methylation-display` | `binned` | Heatmap methylation matrix: `binned` for 50 bp bins or `smoothed` for base-resolution smoothed values. |
| `--signal-colors` | automatic | Optional metaplot colors as `SAMPLE:SIGNAL=COLOR` entries. Signal can be `nuc`, `nucleosome`, `met`, or `5mC`, for example `WT:nuc=#63B8FF WT:5mC=lightcoral KO:nuc=navy KO:5mC=darkred`. |

For CGI, the centered workflow uses a fixed +/-2.5 kb window around the annotation
midpoint. `Promoter`, `Bivalent_Promoter`, and `PRC_Promoter` use a
strand-aware TSS center when a matching 9-column BED is provided. Gene-body
classes use midpoint centering but still use the BED strand for minus-strand
flipping. Without `--bed`, strand is treated as unknown, so promoter plots fall
back to the plus/unknown TSS convention.

The centered metaplot supports any number of samples. If `--signal-colors` is not
provided, PACKAGE assigns a blue gradient to nucleosome profiles and a red
gradient to 5mC profiles in the order supplied to `--samples`.
