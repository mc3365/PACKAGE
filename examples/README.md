# Examples

Runnable examples demonstrating common workflows.

- `ont_region_plot.py`: query an ONT HDF5 database and save a single-molecule plot
- `ont_feature_ecdf.py`: export global per-fiber feature fractions and save ECDF plots
- `ont_annotation_heatmap.py`: export annotation-centered 5mC/nucleosome matrices
  and save heatmap/metaplot figures

Install visualization dependencies before running:

```bash
pip install -e ".[viz]"
python examples/ont_region_plot.py --help
python examples/ont_feature_ecdf.py --help
python examples/ont_annotation_heatmap.py --help
```
