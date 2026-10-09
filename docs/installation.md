# Installation

## Requirements

- Python ≥ 3.10
- Linux or macOS (Windows untested)
- An aligned, coordinate-sorted BAM compatible with the selected workflow
- ONT: MM/ML modification tags and a matching reference FASTA
- PacBio: a fibertools-compatible Fiber-seq BAM; use a FIRE-annotated BAM when
  extracting the optional co-accessibility layer

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

This installs the versions used for the validated ONT workflow: samtools 1.22.1,
modkit 0.5.0, and fibertools-rs 0.8.0. PacBio extraction also uses fibertools-rs.
The `[viz]` extra installs Matplotlib.

## Verify the install

```bash
PACKAGE --version
pytest
```

Both commands should succeed.

## On an HPC cluster

If your cluster already provides these tools in a separate environment, install the
Python package in your analysis environment and put the long-read tool binaries on
`PATH`:

```bash
export PATH=/path/to/ont_env/bin:$PATH
pip install -e ".[dev,viz]"
```

Absolute executable paths can also be set under `extraction` in the YAML config.
