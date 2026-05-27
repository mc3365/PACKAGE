# SLURM scripts

HPC job submission scripts for running `PACKAGE` on the Yale McCleary/Bouchet cluster.

These are **not installed by pip** — they live here as reference and for direct submission
via `sbatch`. They wrap the `PACKAGE` CLI with cluster-specific paths, resource
requests, and module loading.

When porting from the legacy V8 scripts in your scratch directory:
1. Keep the SLURM header (`#SBATCH` directives) as-is
2. Replace the inline `python build_fiber_database_v8.py ...` call with:
   ```bash
   PACKAGE build --config $CONFIG_DIR/your_config.yaml --platform ont
   ```
3. Replace any direct imports of the old scripts with CLI calls
