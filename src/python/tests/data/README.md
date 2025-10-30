# Test Data Directory

This directory contains sample data and mock responses for testing the PySCF Front backend.

## Files

- `mock_pubchem_response.json`: Mock response data from PubChem API
- `sample_h2.xyz`: Sample hydrogen molecule coordinates (for simple, fast calculations)
- `sample_water.xyz`: Sample water molecule coordinates (for testing polyatomic molecules)

## Usage

These files can be loaded in tests using fixtures or standard file operations:

```python
import json
from pathlib import Path

# Load mock API response
data_dir = Path(__file__).parent / 'data'
with open(data_dir / 'mock_pubchem_response.json') as f:
    mock_data = json.load(f)

# Load sample XYZ file
with open(data_dir / 'sample_h2.xyz') as f:
    h2_xyz = f.read()
```

## Adding New Test Data

When adding new test data:
1. Use descriptive filenames
2. Include a comment in this README
3. Keep file sizes small (< 1 KB when possible)
4. Prefer realistic but minimal examples
