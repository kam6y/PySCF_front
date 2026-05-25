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

echo "Cleaning electron-vite and packaging artifacts..."
clean_directory_contents dist
clean_directory_contents release

echo "Cleaning development artifacts from src/python..."

find src/python -type d -name "__pycache__" -exec rm -rf {} + 2>/dev/null || true
find src/python -type f \( -name "*.pyc" -o -name "*.pyo" -o -name "*.pyd" \) -delete 2>/dev/null || true
find src/python -type d -name ".pytest_cache" -exec rm -rf {} + 2>/dev/null || true
find src/python -type f -name ".coverage*" -delete 2>/dev/null || true
find src/python -type d -name "*.egg-info" -exec rm -rf {} + 2>/dev/null || true
find src/python -type d -name ".tox" -exec rm -rf {} + 2>/dev/null || true

echo "✓ Cleanup completed"
