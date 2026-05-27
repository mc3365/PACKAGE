# Installation

## Requirements

- Python ≥ 3.10
- Linux or macOS (Windows untested)

## From source (recommended during development)

```bash
git clone https://github.com/YOUR_USERNAME/PACKAGE.git
cd PACKAGE
pip install -e ".[dev]"
```

The `-e` flag is for *editable* install: changes to the source take effect without
reinstalling.

The `[dev]` extras include pytest, ruff, and mypy.

## Verify the install

```bash
PACKAGE --version
pytest
```

Both commands should succeed.

## On an HPC cluster

If you don't have admin rights, use a virtualenv or conda environment:

```bash
python -m venv ~/PACKAGE_env
source ~/PACKAGE_env/bin/activate
pip install -e ".[dev]"
```
