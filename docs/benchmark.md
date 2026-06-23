# ONT Benchmark

The ONT benchmark records database scale and measures regional-query performance on
an existing PACKAGE HDF5 database. It does not rebuild the database; full-build timing
is recorded separately through the Slurm workflow.

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
