#!/bin/bash
set -e

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
python -c "import fastapi; print(f\"FastAPI: {fastapi.__version__}\")"
python -c "import uvicorn; print(f\"Uvicorn: {uvicorn.__version__}\")"
python -c "import socketio; print('python-socketio: available')"
python -c "import pyscf; print(f\"PySCF: {pyscf.__version__}\")"

echo "=== Verifying packaged FastAPI backend source ==="
python -c "import app; assert hasattr(app, 'app'); print('FastAPI ASGI app import successful')"

echo "✓ Python backend source verification for Linux completed successfully"
