#!/bin/bash
set -e

clean_directory_contents() {
  local target="$1"

  if [ -d "$target" ]; then
    find "$target" -mindepth 1 -maxdepth 1 -exec rm -rf {} +
  else
    mkdir -p "$target"
  fi
}

# Cleanup previous build without deleting a possible Docker mount point
clean_directory_contents conda_env

# Detect Conda base
CONDA_BASE=$(conda info --base 2>/dev/null || echo "$HOME/miniforge3")
source "$CONDA_BASE/etc/profile.d/conda.sh"

# Activate environment
conda activate pyscf-env

# Pack environment
echo "Packing conda environment..."
conda-pack -o conda_env.tar.gz

# Unpack and cleanup
mkdir -p conda_env
tar -xzf conda_env.tar.gz -C conda_env
rm conda_env.tar.gz

echo "Deferring conda-unpack until the packaged app copies conda_env to its runtime location..."

# Clean up development files to reduce size
echo "Cleaning up development files from conda_env..."
find conda_env -type d -name "__pycache__" -exec rm -rf {} + 2>/dev/null || true
find conda_env -type f \( -name "*.pyc" -o -name "*.pyo" -o -name "*.pyd" \) -delete
find conda_env/lib/python*/site-packages -maxdepth 2 -type d \( -name "tests" -o -name "test" \) -exec rm -rf {} + 2>/dev/null || true
find conda_env -type f -name ".coverage*" -delete 2>/dev/null || true
find conda_env -name "*.egg-info" -type d -exec rm -rf {} + 2>/dev/null || true

echo "✓ Conda environment packed and cleaned successfully"
