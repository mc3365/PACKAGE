# Tutorial

PACKAGE uses the same HDF5 schema for ONT and PacBio after platform-specific
extraction outputs have been converted into common intermediate files.

## ONT Workflow

### Prepare the configuration

Copy `configs/ont_template.yaml`, replace every `/path/to` value, and keep the 5mC
and 5hmC layer paths identical. The shared modkit table contains both call types.
The configured methylation threshold is `0.5`; binarization happens during the HDF5
build while raw probabilities are retained.

The BAM must be coordinate sorted, indexed, aligned to the configured reference, and
contain MM/ML tags. PACKAGE checks these requirements before extraction.

For all YAML fields and command-line options, see the
[Parameter Reference](parameters.md).

### Extract and build

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

### Query a region

```bash
PACKAGE query \
  --db /path/to/output/fiber_database.h5 \
  --region chr1:3000000-5000000 \
  --sample sample1
```

## PacBio Phase 1 Workflow

PacBio support currently covers the shared Phase 1 layers:

- nucleosomes;
- 5mC;
- 6mA;
- MSPs.

PacBio 5hmC is not produced by this pathway. FIRE/co-accessibility scores are
validated from `ft extract --all` when present, but they are not yet stored as a
dedicated HDF5 layer.

### Prepare raw fibertools outputs

If starting from a PacBio Fiber-seq or FIRE BAM, run fibertools extraction. A
typical command is:

```bash
ft extract \
  --nuc pacbio_nuc.bed.gz \
  --msp pacbio_msp.bed.gz \
  --m6a pacbio_6ma.bed.gz \
  --cpg pacbio_5mc.bed.gz \
  --all pacbio_all.tsv.gz \
  --simplify \
  yaleFiberAug19_2025.fire.bam
```

`pacbio_5mc.bed.gz` is useful for inspection, but PACKAGE uses
`pacbio_all.tsv.gz` for the packaged 5mC layer because that table contains named
`ref_5mC` and `5mC_qual` fields. The PacBio normalizer skips missing reference
positions such as `-1` and `.`.

### Normalize PacBio outputs for PACKAGE build

The normalization step is file-format normalization, not biological signal
normalization. It checks the raw fibertools files and writes the same
intermediate formats used by the HDF5 builder:

| Raw PacBio file | Builder-ready file |
| --- | --- |
| `pacbio_all.tsv.gz` | modkit-like 5mC TSV/TSV.GZ |
| `pacbio_6ma.bed.gz` | BED12-style 6mA |
| `pacbio_msp.bed.gz` | BED12-style MSP |
| `pacbio_nuc.bed.gz` | flattened nucleosome CSV |

For an already extracted test directory:

```bash
PACBIO_DIR=/path/to/pacbio_mm
OUT=$PACBIO_DIR/package_test
mkdir -p "$OUT"

python - <<PY
from pathlib import Path
from PACKAGE.extract.pacbio import (
    inspect_pacbio_all,
    convert_pacbio_all_5mc_to_modkit,
    normalize_fibertools_block_bed,
)
from PACKAGE.extract.ont import flatten_nucleosome_bed12

pacbio = Path("$PACBIO_DIR")
out = Path("$OUT")

print(inspect_pacbio_all(pacbio / "pacbio_all.tsv.gz"))

print(convert_pacbio_all_5mc_to_modkit(
    pacbio / "pacbio_all.tsv.gz",
    out / "pacbio_5mc_for_PACKAGE.tsv.gz",
))

print(normalize_fibertools_block_bed(
    pacbio / "pacbio_6ma.bed.gz",
    out / "pacbio_6ma_for_PACKAGE.bed",
    expected="single_base",
))

print(normalize_fibertools_block_bed(
    pacbio / "pacbio_msp.bed.gz",
    out / "pacbio_msp_for_PACKAGE.bed",
    expected="interval",
))

normalize_fibertools_block_bed(
    pacbio / "pacbio_nuc.bed.gz",
    out / "pacbio_nuc_for_PACKAGE.bed",
    expected="interval",
)
n = flatten_nucleosome_bed12(
    out / "pacbio_nuc_for_PACKAGE.bed",
    out / "pacbio_nuc_features.csv",
)
print({"nucleosome_rows": n})
PY
```

Check the normalized BED files before building:

```bash
python - <<'PY'
for p in [
    "pacbio_msp_for_PACKAGE.bed",
    "pacbio_6ma_for_PACKAGE.bed",
    "pacbio_nuc_for_PACKAGE.bed",
]:
    print("\nChecking", p)
    with open(p) as f:
        for i, line in zip(range(1, 1001), f):
            c = line.rstrip("\n").split("\t")
            assert len(c) == 12, (p, i, len(c))
            block_count = int(c[9])
            sizes = [x for x in c[10].rstrip(",").split(",") if x]
            starts = [x for x in c[11].rstrip(",").split(",") if x]
            assert block_count == len(sizes) == len(starts), (
                p, i, block_count, len(sizes), len(starts)
            )
    print("PASS")
PY
```

### Build and query a PacBio HDF5

Start from `configs/pacbio_template.yaml`, then point the layer paths to the
normalized files:

```yaml
samples:
  - name: pacbio_test
    layers:
      nucleosomes: /path/to/package_test/pacbio_nuc_features.csv
      5mC: /path/to/package_test/pacbio_5mc_for_PACKAGE.tsv.gz
      6mA: /path/to/package_test/pacbio_6ma_for_PACKAGE.bed
      msp: /path/to/package_test/pacbio_msp_for_PACKAGE.bed

annotations:
  master: /path/to/master_annotations_basic.uniqueID.bed
```

Build and inspect:

```bash
PACKAGE build --config configs/my_pacbio.yaml
PACKAGE info /path/to/output/pacbio_fiber_database.h5
```

Then query a coordinate window:

```bash
PACKAGE query \
  --db /path/to/output/pacbio_fiber_database.h5 \
  --sample pacbio_test \
  --region chr1:3000000-3050000
```

And test per-fiber accessors:

```bash
python - <<'PY'
from PACKAGE import FiberDatabase

db = "/path/to/output/pacbio_fiber_database.h5"
sample = "pacbio_test"
chrom = "chr1"

with FiberDatabase(db) as fdb:
    fibers = fdb.get_fibers_at(chrom, 3000000, 3050000, sample=sample)
    print("n fibers:", len(fibers))
    first = fibers[0] if fibers else None
    print("first fiber:", first)

    if first:
        print("nucleosomes:", fdb.get_nucleosomes(first, chrom, sample=sample))
        print("5mC:", fdb.get_methylation(first, chrom, mod_type="5mC", sample=sample))
        print("6mA:", fdb.get_methylation(first, chrom, mod_type="6mA", sample=sample))
        print("MSP:", fdb.get_msp(first, chrom, sample=sample))
PY
```

### Annotation queries

Annotation BED files are written during `PACKAGE build`. If the HDF5 was built
without annotations, rebuild with the `annotations.master` path added to the YAML.
The package does not yet expose a safe command for adding annotations to an
existing HDF5 file in place.

After rebuilding with annotations:

```bash
python - <<'PY'
from PACKAGE import FiberDatabase

db = "/path/to/output/pacbio_fiber_database.h5"

with FiberDatabase(db) as fdb:
    print(fdb.list_annotations())
    result = fdb.query_annotation_fast(
        "CGI",
        sample="pacbio_test",
        max_regions=50,
        layers=["nucleosomes", "5mC", "6mA", "msp"],
    )
    print(result.shape)
    print(result.head())
PY
```

### Current update behavior

HDF5 technically supports append-mode updates, but PACKAGE currently treats a
database build as a reproducible artifact: annotations and molecular layers are
written during `PACKAGE build`. To add annotations, rebuild from the same
intermediate files with `annotations.master` included.

The same applies to future FIRE/co-accessibility storage. If a new FIRE layer is
added to the schema, existing HDF5 files will need either a rebuild or a dedicated
migration/append command. That append command is not implemented yet because it
needs careful validation to avoid silently mixing layers from different extraction
runs.

## How a Query Uses the Database

For annotation-centered queries such as "find CGI methylation in one sample",
PACKAGE first loads the shared annotation regions, groups them by chromosome, and
then works chromosome by chromosome. The chromosome-level arrays are loaded once
and reused across all regions on that chromosome.

Within a chromosome, PACKAGE finds overlapping fibers from the fiber metadata
arrays, then uses `_indices` to jump directly to each fiber's row range in the
requested feature arrays. For example, the 5mC slice index maps a fiber integer
ID to the start and end rows for that fiber's CpG calls. PACKAGE slices only that
range and masks it to the query interval, then returns per-fiber summaries such
as CpG count and percent methylated.

![PACKAGE CGI methylation query walkthrough](figures/architecture/package_query_walkthrough_cgi_methylation.png)

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

For publication-style plots, provide explicit sample colors:

```bash
python examples/ont_feature_ecdf.py \
  --db /path/to/output/fiber_database.h5 \
  --samples WT KO \
  --sample-colors 'WT=#3b6fb6' 'KO=#c74f46' \
  --outdir figures/ecdf_WT_KO
```

### Centered Heatmap and Metaplot

For fixed-window center/TSS plots, use the centered workflow. Pass the matching
9-column annotation BED when promoter or gene-body strand orientation should be
respected. The BED should be the same annotation file used to build the database,
because strand metadata is matched by chromosome, start, end, and unique region
ID.

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
annotation midpoint. For `Promoter`, `Bivalent_Promoter`, and `PRC_Promoter`,
the window is centered at the strand-aware TSS when the BED contains strand
information. Gene-body classes use a midpoint center but still use strand
information for minus-strand flipping. Minus-strand directional regions are
flipped before averaging so upstream remains on the left.

The heatmap rows are retained annotation regions after averaging across fibers
that span each region. The left heatmap shows nucleosome occupancy, the right
heatmap shows 5mC, and the side bars show row means. With
`--methylation-display binned`, 5mC is shown in 50 bp bins, which gives a smoother
matrix for sparse CpG calls.

The metaplot accepts any number of samples. By default, nucleosome lines use a
blue series and 5mC lines use a red series in the order supplied to `--samples`.
For exact figure colors, override individual sample/signal pairs:

```bash
python examples/ont_centered_heatmap.py \
  --db /path/to/output/fiber_database.h5 \
  --bed /path/to/master_annotations_v4.uniqueID.bed \
  --annotation Promoter \
  --samples lif_d0 lif_d4 \
  --signal-colors 'lif_d0:nuc=#63B8FF' 'lif_d0:5mC=lightcoral' 'lif_d4:nuc=navy' 'lif_d4:5mC=darkred' \
  --methylation-display binned \
  --outdir figures/Promoter_center
```

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

The figure shows each fiber as one row. The thin black line marks the displayed
span of the fiber, gray blocks mark nucleosomes, teal outlines mark MSPs, and
vertical ticks mark the selected modification layer. For 5mC and 5hmC, red ticks
are modified calls and pale gray ticks are unmodified calls.

Use track modes to make the plot less crowded:

- `modification`: selected modification layer only;
- `chromatin`: nucleosomes and MSPs only;
- `full`: nucleosomes, MSPs, and the selected modification layer.

First, a modification-only view is useful when the main question is the
single-molecule distribution of 5mC or 5hmC along each fiber:

```bash
python examples/ont_region_plot.py \
  --db /path/to/output/fiber_database.h5 \
  --region chr1:3000000-3050000 \
  --sample d0 \
  --layer 5mC \
  --tracks modification \
  --max-fibers 100 \
  --out figures/single_molecule/d0_chr1_5mC_only_with_fiber_edges.png
```

![Single-molecule 5mC-only view](figures/visualization/single_molecule/d0_chr1_5mC_only_with_fiber_edges.png)

A chromatin-only view is useful when nucleosome and MSP patterns become hard to
see under base-level modification ticks:

```bash
python examples/ont_region_plot.py \
  --db /path/to/output/fiber_database.h5 \
  --region chr1:3000000-3050000 \
  --sample d0 \
  --tracks chromatin \
  --max-fibers 100 \
  --out figures/single_molecule/d0_chr1_chromatin_only.png
```

![Single-molecule chromatin-only view](figures/visualization/single_molecule/d0_chr1_chromatin_only.png)

A full overlay is useful for inspecting whether methylation, MSPs, and
nucleosomes co-occur on the same individual molecules:

```bash
python examples/ont_region_plot.py \
  --db /path/to/output/fiber_database.h5 \
  --region chr1:3000000-3050000 \
  --sample d0 \
  --layer 5mC \
  --tracks full \
  --max-fibers 100 \
  --out figures/single_molecule/d0_chr1_full.png
```

![Single-molecule full view](figures/visualization/single_molecule/d0_chr1_3000000_3050000_5mC_clean.png)

Limit single-molecule plots to focused windows and use `--max-fibers` for
legibility.
