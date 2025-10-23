"""
Science Analyst Tool Module

This module provides tools for the Science Analyst agent to access and analyze
quantum chemistry calculation results for report generation.

The Science Analyst specializes in:
- Retrieving and interpreting calculation results for a specific calculation ID
- Molecular orbital analysis and visualization
- Spectroscopy data analysis (IR spectra)
- Creating comprehensive scientific reports

All tools are imported from the quantum_calculator.tools module to maintain
a single source of truth and avoid code duplication.

Note: The Supervisor agent is responsible for finding and identifying calculations.
The Science Analyst works with specific calculation IDs provided by the Supervisor.
"""

# Import calculation result retrieval tools
from agent.quantum_calculator.tools import (
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
