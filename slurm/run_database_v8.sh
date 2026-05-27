#!/bin/bash
#SBATCH --job-name=build_fiber_v8
#SBATCH --partition=bigmem
#SBATCH --output=build_db_v8_%j.out
#SBATCH --error=build_db_v8_%j.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=900G
#SBATCH --time=1-00:00:00
#SBATCH --mail-user=meiying.cui@yale.edu
#SBATCH --mail-type=ALL

PYTHON=/vast/palmer/pi/zsmith/mc3365/env/python_env/bin/python
CODE_DIR=/vast/palmer/scratch/zsmith/mc3365/long_read/first_run/code/build_db
SCRIPT=${CODE_DIR}/build_fiber_database_v8.py
DB_PATH=/vast/palmer/scratch/zsmith/mc3365/long_read/first_run/out/h5/fiber_database_v8.h5

echo "Python: $PYTHON"
$PYTHON --version
$PYTHON -c "import h5py, pandas, numpy, tqdm; print('All packages found')"

echo "======================================"
echo "V8 unified database build (integer IDs)"
echo "======================================"

$PYTHON -u $SCRIPT --samples d0 d4

echo "======================================"
echo "Building spatial index..."
echo "======================================"

cd $CODE_DIR
$PYTHON -u -c "
import sys; sys.path.insert(0, '.')
from fiber_database_v8 import FiberDatabase
db = FiberDatabase('${DB_PATH}', build_index=True)
db.close()
"

echo "======================================"
echo "Running validation..."
echo "======================================"

$PYTHON -u test_fiber_database_v8.py --db $DB_PATH

echo "All done!"