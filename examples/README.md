# Examples

Runnable examples demonstrating common workflows.

- `ont_region_plot.py`: query an ONT HDF5 database and save a single-molecule plot
- `ont_feature_ecdf.py`: export global per-fiber feature fractions and save ECDF plots
- `ont_annotation_heatmap.py`: export annotation-centered 5mC/nucleosome matrices
  and save heatmap/metaplot figures
- `ont_centered_heatmap.py`: export fixed-window center/TSS matrices and save
  paired heatmap/metaplot figures
- `ont_benchmark_summary.py`: turn storage, random-query, and annotation-query
  benchmark tables into summary figures

PacBio FIRE co-accessibility uses the main CLI rather than a separate example script:

```bash
mei-fiber coaccess cov --help
mei-fiber coaccess object --help
mei-fiber coaccess rank --help
python -m mei_fiber.benchmark.coaccessibility --help
```

ECDF plots accept `--sample-colors SAMPLE=COLOR`; centered metaplots accept
`--signal-colors SAMPLE:SIGNAL=COLOR` for exact multi-sample figure palettes.

Install visualization dependencies before running:

```bash
pip install -e ".[viz]"
python examples/ont_region_plot.py --help
python examples/ont_feature_ecdf.py --help
python examples/ont_annotation_heatmap.py --help
python examples/ont_centered_heatmap.py --help
```
