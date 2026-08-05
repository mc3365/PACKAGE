# Step-by-step ONT workflow

This guide starts with an aligned Oxford Nanopore BAM and ends with a queryable
PACKAGE database and a single-molecule regional plot. It follows the workflow tested
with the d0 and d4 mouse embryonic stem-cell data.

This guide focuses on the ONT workflow. PacBio Phase 1 support is documented in
[`docs/tutorial.md`](docs/tutorial.md) and [`configs/pacbio_template.yaml`](configs/pacbio_template.yaml).

For a complete list of command-line and YAML options, see
[`docs/parameters.md`](docs/parameters.md).

## 1. Clone the repository

Choose a directory where you keep code projects, then clone PACKAGE:

```bash
git clone https://github.com/mc3365/PACKAGE.git
cd PACKAGE
```

To update an existing clone:

```bash
git checkout main
git pull --ff-only
```

## 2. Install PACKAGE and the ONT tools

The supplied Conda environment includes the versions used during validation:

- samtools 1.22.1
- modkit 0.5.0
- fibertools-rs 0.8.0

Create the environment and install PACKAGE:

```bash
conda env create -f environments/ont.yml
conda activate package-ont
pip install -e ".[viz]"
```

Confirm that the command-line tools are available:

```bash
PACKAGE --version
samtools --version | head -1
modkit --version
ft --version
```

On an HPC system, the Python package and ONT tools may live in separate environments.
In that case, set `modkit_executable`, `ft_executable`, and `samtools_executable` in
the YAML configuration to their absolute paths.

## 3. Check the input BAM

The input BAM must:

- be aligned and coordinate sorted;
- have a BAM index (`.bai`);
- use chromosome names compatible with the reference FASTA; and
- contain paired MM and ML modification tags.

PACKAGE checks these requirements before extraction. Basic checks can also be run
directly:

```bash
samtools quickcheck sample.bam
samtools view -H sample.bam | grep '^@HD'
samtools view sample.bam | head -1
```

The BAM and reference FASTA must describe the same genome assembly. For the validated
mouse workflow, the reference path was an mm10 FASTA.

## 4. Prepare the configuration

Make a project-specific copy of the ONT template:

```bash
cp configs/ont_template.yaml configs/my_ont.yaml
```

Edit these fields in `configs/my_ont.yaml`:

- `output_dir`: directory for the HDF5 database;
- `output_file`: HDF5 filename;
- `reference`: reference FASTA used for alignment;
- `samples[].name`: short, unique sample name;
- `samples[].bam`: aligned, modification-tagged BAM;
- `samples[].layers`: paths for the extraction outputs; and
- `extraction`: threads, overwrite behavior, and executable paths.

The 5mC and 5hmC entries must point to the same modkit TSV because both call types are
stored in that file. Keep `methylation_threshold: 0.5` to reproduce the validated
binarization rule. Raw probabilities are also retained in the HDF5 database.

For a second sample, add another item under `samples` with its own BAM and output
paths. Sample names become top-level groups in the database.

Common YAML settings:

| Field | Typical value | When to change it |
| --- | --- | --- |
| `parameters.methylation_threshold` | `0.5` | Keep this value to reproduce the validated ONT binarization. |
| `build.build_spatial_index` | `true` | Keep enabled for fast coordinate queries. |
| `extraction.threads` | `16` | Match available CPU cores. |
| `extraction.overwrite` | `false` | Set to `true` only when regenerating extraction products. |
| `extraction.validate_bam` | `true` | Keep enabled unless doing a very controlled rerun. |
| `extraction.write_manifest` | `true` | Keep enabled for provenance and troubleshooting. |
| `extraction.modkit_executable` | `modkit` or absolute path | Use an absolute path on HPC if tools live outside the Python environment. |
| `extraction.ft_executable` | `ft` or absolute path | Use an absolute path on HPC if needed. |
| `extraction.samtools_executable` | `samtools` or absolute path | Use an absolute path on HPC if needed. |

## 5. Run a small extraction test

Before processing a full BAM, it is useful to test a small genomic interval. This
checks the environment, reference, BAM tags, and output paths at low cost:

```bash
samtools view -@ 4 -h -b sample.bam chr1:3000000-5000000 \
  | samtools sort -@ 4 -o sample_test.bam
samtools index sample_test.bam
samtools view -c sample_test.bam
```

Create a second YAML file that points to `sample_test.bam` and writes to a separate
test directory. Then run:

```bash
PACKAGE extract \
  --platform ont \
  --config configs/ont_smoke_test.yaml \
  --samples sample_test
```

Extraction command options:

| Option | Meaning |
| --- | --- |
| `--platform ont` | Use the ONT extractor. |
| `--config` | YAML configuration file. |
| `--samples` | Optional sample filter. Repeat the flag to process multiple selected samples. |

Successful extraction produces:

- a shared modkit TSV for 5mC and 5hmC;
- BED output for 6mA;
- BED output for MSPs;
- a flattened nucleosome CSV; and
- `PACKAGE_manifest_<sample>.json` with BAM QC, tool versions, commands, and output
  sizes.

Inspect the manifest before continuing:

```bash
python -m json.tool /path/to/PACKAGE_manifest_sample_test.json | less
```

## 6. Extract the full samples

Run full extraction on a compute node or submit it as a Slurm job. To process every
sample in the YAML:

```bash
PACKAGE extract --platform ont --config configs/my_ont.yaml
```

To process selected samples, repeat `--samples`:

```bash
PACKAGE extract \
  --platform ont \
  --config configs/my_ont.yaml \
  --samples d0 \
  --samples d4
```

With `overwrite: false`, valid existing outputs are reused. Set `overwrite: true` only
when the extraction products should be regenerated.

## 7. Build the HDF5 database

Run the unit tests before a full build when developing from a Git checkout:

```bash
python -m pytest tests/test_utils.py tests/test_builder.py -q
```

Build every configured sample:

```bash
PACKAGE build --config configs/my_ont.yaml
```

Build command options:

| Option | Meaning |
| --- | --- |
| `--config` | YAML file describing extracted layer paths and database output. |
| `--samples` | Optional sample filter. Repeat to build multiple selected samples. |
| `--build-index` | Build the spatial-index sidecar after HDF5 creation. |
| `--no-build-index` | Skip spatial-index creation for this run. |

The builder creates `output_dir/output_file`. When `build.build_spatial_index: true`
is set in the YAML, it also creates a cached spatial-index sidecar named
`<output>.index.pkl` for fast region queries. Keep that file beside the HDF5 database.
The builder opens the HDF5 path for writing, so use a new filename or preserve the
previous HDF5 file before rebuilding.

Inspect the completed database:

```bash
PACKAGE info /path/to/output/fiber_database.h5
ls -lh /path/to/output/fiber_database.h5
ls -lh /path/to/output/fiber_database.index.pkl
```

Confirm that the expected samples and annotations appear and that the fiber counts are
consistent with the extraction inputs.

## 8. Query a region

Use the command line for a quick overlap check:

```bash
PACKAGE query \
  --db /path/to/output/fiber_database.h5 \
  --region chr1:3000000-5000000 \
  --sample d0
```

Query command options:

| Option | Meaning |
| --- | --- |
| `--db` | HDF5 database path. |
| `--region` | Region in `chrom:start-end` form. |
| `--sample` | Optional sample filter. |
| `--out` | Optional text file for overlapping fiber IDs. |

For repeated regional analysis, create the spatial index once:

```python
from PACKAGE.db import FiberDatabase

with FiberDatabase("/path/to/output/fiber_database.h5", build_index=True) as db:
    print(db.get_summary())
```

The index is saved beside the HDF5 database and loaded by later sessions.

## 9. Plot individual molecules

The example plot shows nucleosomes, MSPs, and one modification layer across fibers
overlapping a region:

```bash
python examples/ont_region_plot.py \
  --db /path/to/output/fiber_database.h5 \
  --region chr1:3000000-3050000 \
  --sample d0 \
  --layer 5mC \
  --tracks full \
  --max-fibers 100 \
  --hide-unmodified \
  --out d0_chr1_region.png
```

Use a focused interval, usually 10-50 kb, so individual molecular patterns remain
readable. The same command accepts `5hmC` or `6mA` as the modification layer.

Region plot options:

| Option | Meaning |
| --- | --- |
| `--db` | HDF5 database path. |
| `--region` | Focused interval in `chrom:start-end` form. |
| `--sample` | Optional sample name. |
| `--layer` | Modification layer: `5mC`, `5hmC`, or `6mA`. |
| `--tracks` | `full` shows nucleosome/MSP plus the selected modification; `modification` shows only the selected modification; `chromatin` shows nucleosome/MSP only. |
| `--max-fibers` | Maximum number of fiber rows shown. |
| `--hide-unmodified` | Hide unmodified 5mC/5hmC calls for a cleaner modification-focused view. |
| `--out` | Output figure path. |

## 10. Plot summary figures

Global ECDF plots summarize per-fiber feature fractions:

```bash
python examples/ont_feature_ecdf.py \
  --db /path/to/output/fiber_database.h5 \
  --samples d0 d4 \
  --outdir figures/ecdf
```

Useful options:

| Option | Meaning |
| --- | --- |
| `--samples` | One or more samples. Omit to include all samples. |
| `--max-fibers-per-chrom` | Optional downsampling cap for exploratory plots. |

For annotation-body heatmaps and metaplots:

```bash
python examples/ont_annotation_heatmap.py \
  --db /path/to/output/fiber_database.h5 \
  --annotation CGI \
  --sample d0 \
  --max-regions 50 \
  --min-fibers 10 \
  --outdir figures/CGI_d0_test
```

Omit `--max-regions` for full annotation plots. Keep it for smoke tests.

For fixed-window center/TSS heatmaps and metaplots:

```bash
python examples/ont_centered_heatmap.py \
  --db /path/to/output/fiber_database.h5 \
  --bed /path/to/master_annotations_v4.uniqueID.bed \
  --annotation CGI \
  --samples d0 d4 \
  --methylation-display binned \
  --min-fibers 10 \
  --outdir figures/CGI_center
```

Centered plot options:

| Option | Meaning |
| --- | --- |
| `--annotation` | Annotation class, such as `CGI`. |
| `--samples` | One or more samples to plot together in the metaplot. |
| `--bed` | Optional 9-column annotation BED for strand-aware TSS/gene orientation. |
| `--max-regions` | Optional smoke-test limit. Omit for full figures. |
| `--min-fibers` | Minimum spanning fibers required for a region. Default is `10`. |
| `--methylation-display` | `binned` for 50 bp bins or `smoothed` for smoothed base-resolution 5mC. |

For more detail on plotting outputs and defaults, see
[`docs/parameters.md`](docs/parameters.md#plotting).

## 11. Run benchmarks

Use the annotation benchmark for a biologically meaningful speed test:

```bash
python -m PACKAGE.benchmark.annotation \
  --db /path/to/output/fiber_database.h5 \
  --sample d0 \
  --annotation CGI \
  --layers nucleosomes 5mC 5hmC msp \
  --max-regions 50 \
  --repeats 3 \
  --warmups 1 \
  --outdir benchmark/d0_CGI
```

Use the random-region benchmark for technical scaling across window sizes:

```bash
python -m PACKAGE.benchmark.ont \
  --db /path/to/output/fiber_database.h5 \
  --config configs/my_ont.yaml \
  --sample d0 \
  --outdir benchmark/d0 \
  --region-sizes-kb 1 10 100 1000 \
  --n-regions 25 \
  --repeats 3 \
  --query-modes spatial_index \
  --build-index
```

For full production databases, avoid the slow array-scan baseline unless you
specifically need it. Use `--query-modes spatial_index` for routine benchmarking.

Benchmark options and outputs are described in
[`docs/benchmark.md`](docs/benchmark.md) and
[`docs/parameters.md`](docs/parameters.md#random-region-benchmark).

## 12. Verify a development checkout

Contributors should install the development dependencies and run the complete checks:

```bash
pip install -e ".[dev,viz]"
python -m pytest -q
python -m ruff check src tests examples
```

GitHub Actions runs the automated tests after each push. The local checks should pass
before changes are committed.

## Troubleshooting

### `conda: command not found`

Initialize Conda for the current shell or use the full path to the environment's
executables. On batch systems, shell startup files may not be loaded by Slurm.

### An executable is not found

Activate the ONT environment, add its `bin` directory to `PATH`, or set the three
executable paths explicitly in the YAML configuration.

### `gzip.BadGzipFile`

Some older extraction outputs have a `.gz` suffix but contain plain text. PACKAGE can
read both forms, but new modkit outputs should be checked with:

```bash
file sample_raw_mods.tsv.gz
gzip -t sample_raw_mods.tsv.gz
```

### Region queries report that no spatial index was found

The database is still usable, but queries fall back to an array scan. Open the database
once with `build_index=True` as shown in Step 8.

### GitHub authentication fails

Refresh the GitHub CLI login, then retry the push:

```bash
gh auth login -h github.com
git push
```

## Development history

PACKAGE began as a Python package skeleton and was filled in against the most recent
legacy V8 scripts. The database reader was first checked for output parity, followed by
the builder and then the ONT extraction workflow. Real d0 and d4 data were used to
validate the HDF5 build, regional queries, modification binarization, BAM QC, manifests,
and extraction outputs. The validated ONT milestone is tagged `v0.4.0-ont`.
