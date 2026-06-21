# Installation

## Requirements

- Python ≥ 3.10
- Linux or macOS (Windows untested)
- An aligned, coordinate-sorted ONT BAM containing MM/ML modification tags

## From source (recommended during development)

```bash
git clone https://github.com/mc3365/PACKAGE.git
cd PACKAGE
conda env create -f environments/ont.yml
conda activate package-ont
pip install -e ".[dev,viz]"
```

The `-e` flag is for *editable* install: changes to the source take effect without
reinstalling.

This installs the versions validated on HPC: samtools 1.22.1, modkit 0.5.0, and
fibertools-rs 0.8.0. The `[viz]` extra installs Matplotlib.

## Verify the install

```bash
PACKAGE --version
pytest
```

Both commands should succeed.

## On an HPC cluster

If your cluster already provides these tools in a separate environment, install the
Python package in your analysis environment and put the ONT tool binaries on `PATH`:

```bash
export PATH=/path/to/ont_env/bin:$PATH
pip install -e ".[dev,viz]"
```

Absolute executable paths can also be set under `extraction` in the YAML config.
