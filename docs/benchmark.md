# Benchmarking and Validation

The ONT benchmark records database scale and measures regional-query performance on
an existing PACKAGE HDF5 database. It does not rebuild the database; full-build timing
is recorded separately through the Slurm workflow.

PacBio FIRE validation compares legacy and PACKAGE `Cov.bed`, object, pair-ranking,
and stitched-region outputs. These are complementary questions: the ONT benchmarks
measure storage and query performance, while the PacBio validation measures workflow
equivalence and explains intentional differences from the legacy implementation.

## Run the benchmark

Use a compute node and write each sample to a separate output directory. Build the
spatial index during the first run if the database does not already have one:

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

For another sample, reuse the saved index and change the output directory:

```bash
python -m PACKAGE.benchmark.ont \
  --db /path/to/fiber_database.h5 \
  --config configs/my_ont.yaml \
  --sample d4 \
  --outdir benchmark/d4
```

The default design measures 25 fixed-seed, fiber-centered regions at each of four
window sizes, with three timed repeats. Every selected region receives one untimed
warm-up query. This avoids mixing one-time HDF5 and fiber-ID initialization costs into
the regional-query measurements.

For full production databases, use `--query-modes spatial_index`. The array-scan mode
is useful as a correctness and baseline comparison on small databases, but it can take
many hours when repeated across millions of fibers.

## Outputs

- `storage_summary.json`: HDF5 size, spatial-index size, annotations, and per-sample
  record totals;
- `database_record_counts.csv`: per-chromosome fiber and layer counts;
- `intermediate_file_sizes.csv`: unique extraction-product sizes from the YAML;
- `query_speed_results.csv`: every region, repeat, fiber count, and elapsed time;
- `query_speed_summary.csv`: median, quartiles, and p95 by mode and window size;
- `query_speed_metadata.json`: database path, package version, random seed, machine,
  and benchmark design;
- `query_speed.png`: raster preview; and
- `query_speed.pdf`: vector figure for the manuscript.

By default, when a spatial index is available, the same regions are measured with both
the array scan and spatial index. Use `--query-modes spatial_index` or
`--query-modes array_scan` to select a subset.

## Reporting

Report the median and interquartile range rather than a single fastest run. Keep the
raw CSV and metadata JSON with the manuscript analysis so the panel can be regenerated.
The benchmark measures warm-cache regional queries; it does not claim to measure
first-access disk latency.

The random-window benchmark measures coordinate lookup only. Each point in the
random-query scatter is one call to `FiberDatabase.get_fibers_at(chrom, start, end)`,
which returns the fiber IDs overlapping that interval. It does not load methylation,
nucleosome, 6mA, or MSP rows. The point color shows how many fibers were returned.
Very high-overlap regions can be slow because returning hundreds of thousands of
fiber IDs is itself a large output.

## Summary Figures

After running the storage, random-window, and annotation benchmarks, combine the
outputs into manuscript-style summary panels:

```bash
python examples/ont_benchmark_summary.py \
  --storage-summary benchmark/d0/storage_summary.json \
  --query-results benchmark/d0/query_speed_results.csv \
  --annotation-summary benchmark/annotation_d0/annotation_query_summary.csv \
  --bam /path/to/d0.output.fiberseq.bam \
  --bam /path/to/d4.output.fiberseq.bam \
  --outdir figures/benchmark_summary
```

The BAM paths are optional. When supplied, the storage panel shows the aligned BAM
footprint alongside extracted intermediates, the PACKAGE HDF5 file, and the spatial
index. The `intermediate_size_bytes` value in `storage_summary.json` is the sum of
the extracted layer files listed in the YAML; it does not include the original BAM.

The summary script writes:

- `benchmark_storage_summary.*`: BAM/intermediate/HDF5/index size comparison;
- `benchmark_record_counts.*`: per-sample record counts by layer. For example, a
  5mC value of 0.5 means about 0.5 billion CpG records stored for that sample;
- `benchmark_random_query_scatter.*`: random-window query time with returned-fiber
  counts and outliers marked;
- `query_speed_summary_without_outliers.csv`: random-window query summary excluding
  regions with very large returned-fiber counts;
- `benchmark_annotation_queries.*`: annotation query timing;
- `package_ont_workflow.*`: current ONT workflow schematic; and
- `benchmark_summary.md`: key numeric ratios.

## How to Read the Plots

The benchmark summary plots are meant to answer three different questions:
how much data was integrated, how fast coordinate lookup is, and how fast a real
annotation-centered biological query is.

### `benchmark_storage_summary`

This panel compares the data footprint of the ONT workflow. The extracted
intermediate files are the layer files produced after running `modkit` and
`fibertools-rs` such as raw modification TSVs, 6mA BED files, MSP BED files,
and nucleosome CSV files. They are not the original aligned BAM files. When
`--bam` paths are supplied to `examples/ont_benchmark_summary.py`, the BAM sizes
are shown separately.

The HDF5 bar shows the integrated PACKAGE database size, and the spatial-index
bar shows the sidecar index used for fast coordinate lookup. A large difference
between intermediate files and HDF5 means the package has consolidated many
large per-layer files into a smaller queryable database.

![ONT benchmark storage summary](figures/benchmark/benchmark_storage_summary.png)

### `benchmark_record_counts`

This panel shows how many records were stored for each sample and molecular
layer. For example, a 5mC value near 0.5 means about 0.5 billion CpG-level
records are stored for that sample. The stacked bars summarize database scale:
fibers, 5mC calls, 5hmC calls, 6mA calls, MSP intervals, and nucleosome records.

The 5mC and 5hmC counts can be similar because both are parsed from the same
`modkit extract full` CpG table, while the stored probabilities and binary calls
for each modification code are kept as separate layers.

![ONT benchmark record counts](figures/benchmark/benchmark_record_counts.png)

### `benchmark_random_query_scatter`

This panel measures coordinate lookup speed with random genomic windows. Each
dot is one timed call to:

```python
FiberDatabase.get_fibers_at(chrom, start, end)
```

That call returns the IDs of fibers overlapping the interval. It does not load
the raw methylation, nucleosome, 6mA, or MSP rows. The x-axis is the query window
size, the y-axis is elapsed time, and the color shows how many fibers were
returned. This is useful for testing the spatial index itself.

Occasional slow points usually mean the selected window overlaps an unusually
large number of fibers. In that case the time includes returning a very large
result, not just finding the interval. The summary CSV without outliers is useful
for describing typical lookup speed, while the full CSV should be kept for
transparency.

![ONT random-window query benchmark](figures/benchmark/benchmark_random_query_scatter.png)

### `benchmark_annotation_queries`

This panel measures a more biologically realistic query. For each annotation
class, `--max-regions 50` takes the first 50 stored regions of that class, finds
the fibers overlapping those regions, and computes per-fiber summaries for the
requested layers.

The benchmark times:

```python
FiberDatabase.query_annotation_fast(...)
```

For each `(annotation region, overlapping fiber)` row, PACKAGE summarizes the
requested molecular features, such as nucleosome count and spacing, CpG count
and percent methylated, 5hmC count, and MSP count/width. This is not just a
coordinate lookup and it does not return every raw base-level array; it measures
the common analysis pattern: "for this annotation class, find overlapping fibers
and summarize their molecular features."

The row counts shown in the output are therefore the returned
`(annotation region, fiber)` rows, not the number of annotation regions. For
example, 50 promoter regions can produce more than 1,000 returned rows if many
individual fibers overlap those promoters.

![ONT annotation query benchmark](figures/benchmark/benchmark_annotation_queries.png)

### `package_ont_workflow`

This schematic summarizes the ONT branch of PACKAGE: aligned Fiber-seq
BAM files are processed by `modkit` and `fibertools-rs`, the extracted layers are
packed into an indexed HDF5 database, and downstream users query or visualize
single-molecule molecular features. It is retained as an ONT-specific schematic;
the PacBio/FIRE path is documented separately below.

![PACKAGE ONT workflow](figures/benchmark/package_ont_workflow.png)

## PacBio FIRE Co-Accessibility Validation

Run this comparison after generating PACKAGE `Cov.bed`, object, pair-rank, and
cluster-rank outputs. The legacy pair table contains both `A -> B` and `B -> A`, so
the validator first canonicalizes element order and compares one row per biological
pair.

```bash
python -m PACKAGE.benchmark.coaccessibility \
  --legacy-ce /path/to/legacy/ce_rank.txt \
  --package-ce /path/to/package/ce_rank.txt \
  --legacy-cluster /path/to/legacy/cluster_rank.txt \
  --package-cluster /path/to/package/cluster_rank.txt \
  --legacy-cov /path/to/legacy/Cov.sorted.bed \
  --package-cov /path/to/package/Cov.sorted.bed \
  --cov-difference /path/to/Cov.missing_from_PACKAGE.bed \
  --legacy-object /path/to/legacy/scored_obj.json \
  --package-object /path/to/package/scored_PACKAGE_obj.json \
  --outdir benchmark/coaccess_validation
```

The command writes:

- `validation_summary.json`: pair, cluster, Cov, and object-level metrics;
- `pair_agreement.csv` and `cluster_agreement.csv`: scores and ranks for shared
  records;
- `pair_membership.csv` and `cluster_membership.csv`: shared and method-specific
  records;
- `top_k_overlap.csv`: overlap among the strongest `k` constituent pairs; and
- `coaccess_validation.png/.pdf`: a four-panel validation figure.

The `--threshold 0.10` default is the legacy per-fiber FIRE accessibility cutoff:
scores below 0.10 are treated as accessible. It is not a p-value, methylation
threshold, or `Super` cutoff. Pair and cluster `Super` labels are assigned later from
the elbow of their ranked corrected-score distributions.

### Validated reference result

The current reference dataset produced the following comparison:

| Metric | Result | Interpretation |
| --- | ---: | --- |
| Unique constituent pairs | 74,435 legacy; 74,435 PACKAGE; all shared | Exact biological pair membership |
| Pair score Pearson | 0.990 | Corrected score magnitudes are highly concordant |
| Pair/rank Spearman | 0.997 | Pair ordering is nearly identical |
| Pair `Super` agreement | 99.93% | 52 boundary differences among 74,435 pairs |
| Legacy stitched regions recovered | 15,287/15,287 | No legacy region was lost |
| Additional PACKAGE regions | 601 | All match the legacy `0.5` sentinel-collision diagnostic |
| Cluster score Pearson | 0.990 | Shared stitched-region scores are highly concordant |
| Cluster/rank Spearman | 0.986 | Shared region ordering is strongly concordant |
| Cov rows | 6,739,056 legacy; 6,739,053 PACKAGE | Three FIRE-only fibers were absent from the HDF5 molecule table |
| Cov reconstruction | Exact after adding 3 known rows | No unexplained Cov difference remained |

![PacBio FIRE co-accessibility validation](figures/benchmark/coaccess_validation.png)

Read the four panels as follows:

1. **Constituent-pair scores:** points close to the identity line indicate similar
   distance-corrected scores over several orders of magnitude.
2. **Constituent-pair ranks:** small rank shifts arise because the legacy distance fit
   uses directional duplicates while PACKAGE fits unique unordered pairs.
3. **Top-k pair overlap:** 90% of the top 10 and 95-98% of larger top-k sets are
   shared, showing that the strongest biological results are stable.
4. **Shared stitched regions:** most shared region scores lie near identity. The
   low-score horizontal feature is emphasized by the logarithmic axes and reflects
   small low-end scoring differences. PACKAGE-only regions are not drawn in this
   shared-only panel.

The 601 additional PACKAGE regions are not unexplained false positives. The legacy
cluster code used numeric `0.5` both for a legitimate pair result and as an invalid
pair sentinel, so valid exact-0.5 cases could be removed. PACKAGE uses an unambiguous
missing value; all 601 additional regions satisfy the sentinel-collision diagnostic.

Pair membership is identical, while corrected scores and `Super` calls are highly
concordant but not numerically identical because the distance fitting and legacy
sentinel handling differ.

## Annotation Query Benchmark

Use the annotation benchmark for biologically meaningful workflows such as CGI,
promoter, enhancer, or super-enhancer queries. This times
`FiberDatabase.query_annotation_fast`: annotation class -> overlapping fibers ->
per-fiber layer summaries.

```bash
python -m PACKAGE.benchmark.annotation \
  --db /path/to/fiber_database.h5 \
  --sample d0 \
  --annotation CGI \
  --layers nucleosomes 5mC 5hmC 6mA msp \
  --max-regions 50 \
  --repeats 3 \
  --warmups 1 \
  --outdir benchmark/d0_CGI
```

Repeat `--annotation` to benchmark more than one annotation class in the same run.
Outputs are:

- `annotation_query_results.csv`: one row per timed repeat;
- `annotation_query_summary.csv`: median, quartiles, p95, and result-row counts; and
- `annotation_query_metadata.json`: database, sample, layers, annotations, and run
  settings.

This benchmark is usually more interpretable for figures than random windows because
it measures a real analysis pattern. The random-region benchmark remains useful for
technical scaling across interval sizes.

`--max-regions 50` processes the first 50 regions of each annotation class from the
stored annotation table; it is a cap for reproducible smoke benchmarks, not a random
sample. The reported row count is not the number of annotation regions. It is the
number of `(annotation region, overlapping fiber)` rows returned after finding fibers
that overlap those 50 regions. For each row, PACKAGE computes requested layer
summaries such as nucleosome counts, CpG counts, percent methylated, 5hmC counts, and
MSP counts/widths. It does not return all raw base-level arrays in this benchmark.

## Full-build resources

Full-build timing requires a new output HDF5 file. Copy the working YAML, change
`output_file` to a benchmark-specific filename, and submit the instrumented job:

```bash
sbatch --export=ALL,\
CONFIG=/absolute/path/configs/benchmark_ont.yaml,\
PYTHON=/absolute/path/env/bin/python,\
PACKAGE_ROOT=/absolute/path/PACKAGE \
slurm/benchmark_full_build.sh
```

The job refuses to overwrite an existing database. It writes the raw GNU `time`
report and a JSON report containing wall time, peak RSS, Slurm resources, Git commit,
package version, final HDF5 size, and sample fiber counts. See `slurm/README.md` for
the cluster workflow.
