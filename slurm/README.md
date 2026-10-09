# SLURM scripts

HPC job submission scripts for running `MEI-Fiber` on the Yale McCleary/Bouchet cluster.

These are **not installed by pip** — they live here as reference and for direct submission
via `sbatch`. They wrap the `MEI-Fiber` CLI with cluster-specific paths, resource
requests, and module loading.

## Full-build benchmark

`benchmark_full_build.sh` measures a complete HDF5 build with GNU `time` and writes a
JSON report containing the Git commit, package version, Slurm allocation, elapsed time,
peak memory, output size, and sample fiber counts.

First copy the working YAML and give it a new output filename. The benchmark script
refuses to overwrite an existing HDF5 file by default.

Submit with absolute paths:

```bash
sbatch --export=ALL,\
CONFIG=/absolute/path/configs/benchmark_ont.yaml,\
PYTHON=/absolute/path/env/bin/python,\
MEI_FIBER_ROOT=/absolute/path/MEI-Fiber \
slurm/benchmark_full_build.sh
```

Results are written under `benchmark/build_<job_id>/` unless `RESULTS_DIR` is supplied.
The HDF5 database itself is written to the location specified by the benchmark YAML.

`run_database_v8.sh` is retained only as a record of the legacy V8 build. New MEI-Fiber
builds should use the benchmark script or `mei-fiber build --config <yaml>`.
