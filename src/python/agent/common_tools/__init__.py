"""
Common Tools Package

This package provides a centralized collection of tools for AI agents,
organized by their primary function:

- execution_tools: Tools that modify system state (execution, deletion, configuration)
- analysis_tools: Tools that read and analyze data (read-only operations)

This organization ensures clear separation of concerns and makes it easy
for agents to access only the tools appropriate for their role.
"""

# Import execution tools (state-modifying operations)
from .execution_tools import (
    start_quantum_calculation,
    delete_calculation,
    _execute_confirmed_deletion,
    delete_cube_files,
    convert_smiles_to_xyz,
    validate_xyz_format,
    search_pubchem_by_name,
    update_app_settings,
)

# Import analysis tools (read-only operations)
from .analysis_tools import (
    list_all_calculations,
    find_calculations,
    get_calculation_details,
    get_supported_parameters,
    get_app_settings,
    get_system_resources,
    get_molecular_orbitals,
    generate_orbital_cube,
    list_cube_files,
    generate_ir_spectrum,
)

__all__ = [
    # Execution tools
    'start_quantum_calculation',
    'delete_calculation',
    '_execute_confirmed_deletion',
    'delete_cube_files',
    'convert_smiles_to_xyz',
    'validate_xyz_format',
    'search_pubchem_by_name',
    'update_app_settings',

    # Analysis tools
    'list_all_calculations',
    'find_calculations',
    'get_calculation_details',
    'get_supported_parameters',
    'get_app_settings',
    'get_system_resources',
    'get_molecular_orbitals',
    'generate_orbital_cube',
    'list_cube_files',
    'generate_ir_spectrum',
]
