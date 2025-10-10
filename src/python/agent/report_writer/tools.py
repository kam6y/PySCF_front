"""
Report Writer Tool Module

This module provides tools for the Report Writer agent to access and analyze
quantum chemistry calculation results for report generation.

The Report Writer specializes in:
- Retrieving and interpreting calculation results
- Molecular orbital analysis and visualization
- Spectroscopy data analysis (IR spectra)
- Creating comprehensive scientific reports

All tools are imported from the quantum_calc.tools module to maintain
a single source of truth and avoid code duplication.
"""

# Import calculation result retrieval tools
from agent.quantum_calc.tools import (
    get_calculation_details,
    get_molecular_orbitals,
    generate_orbital_cube,
    list_cube_files,
    delete_cube_files,
    generate_ir_spectrum,
)

__all__ = [
    'get_calculation_details',
    'get_molecular_orbitals',
    'generate_orbital_cube',
    'list_cube_files',
    'delete_cube_files',
    'generate_ir_spectrum',
]
