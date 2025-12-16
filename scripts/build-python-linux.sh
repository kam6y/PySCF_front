#!/bin/bash
set -e

# Cleanup previous build
rm -rf python_dist

cd src/python

# Detect Conda base
CONDA_BASE=$(conda info --base 2>/dev/null || echo "$HOME/miniforge3")
export PATH="$CONDA_BASE/envs/pyscf-env/bin:$PATH"
export CONDA_DEFAULT_ENV=pyscf-env

echo "=== Python Environment Verification ==="
which python
python --version

echo "=== Checking critical dependencies ==="
python -c "import gunicorn; print(f\"Gunicorn: {gunicorn.__version__}\")"
python -c "import watchdog; print(\"Watchdog: Available\")"
python -c "import flask; print(f\"Flask: {flask.__version__}\")"
python -c "import pyscf; print(f\"PySCF: {pyscf.__version__}\")"

echo "=== Dependencies verified, starting PyInstaller ==="
pyinstaller --distpath ../../python_dist --workpath ../../build/pyinstaller --noconfirm pyscf_front_api.spec

echo "✓ Python build for Linux completed successfully"
