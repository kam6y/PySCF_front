"""Calculation method defaults and parameter constraints.

This module defines default parameter values for each quantum chemistry calculation method
and validation constraints for parameters. This serves as the single source of truth for
business logic related to calculation parameters.
"""

from typing import Dict, Any

# Method-specific default values
# Each calculation method has its own recommended default parameters based on
# quantum chemistry best practices and computational efficiency considerations
METHOD_DEFAULTS: Dict[str, Dict[str, Any]] = {
    'DFT': {
        'basis_function': '6-31G(d)',
        'exchange_correlation': 'B3LYP',
        'memory_mb': 2000,
        'optimize_geometry': True
    },
    'HF': {
        'basis_function': '6-31G(d)',
        'memory_mb': 2000,
        'optimize_geometry': True
    },
    'MP2': {
        'basis_function': '6-31G(d)',
        'memory_mb': 3000,
        'optimize_geometry': True
    },
    'CCSD': {
        'basis_function': 'cc-pVDZ',  # Correlation-consistent basis recommended for CCSD
        'memory_mb': 4000,            # Higher memory requirements
        'frozen_core': True,          # Frozen core approximation to reduce cost
        'optimize_geometry': False    # Geometry optimization not available
    },
    'CCSD_T': {
        'basis_function': 'cc-pVDZ',  # Correlation-consistent basis recommended
        'memory_mb': 4000,            # Higher memory requirements
        'frozen_core': True,          # Frozen core approximation to reduce cost
        'optimize_geometry': False    # Geometry optimization not available
    },
    'TDDFT': {
        'basis_function': '6-31G(d)',
        'exchange_correlation': 'B3LYP',
        'memory_mb': 2000,
        'tddft_nstates': 10,          # Default number of excited states
        'tddft_method': 'TDDFT',      # Full TDDFT (vs TDA approximation)
        'tddft_analyze_nto': False,   # NTO analysis off by default
        'optimize_geometry': False    # Geometry optimization not available
    },
    'CASCI': {
        'basis_function': '6-31G(d)',
        'memory_mb': 3000,
        'ncas': 4,                    # Number of active orbitals
        'nelecas': 4,                 # Number of active electrons
        'natorb': True,               # Natural orbital transformation recommended
        'max_cycle_micro': 3,         # CI solver iterations
        'optimize_geometry': False    # Geometry optimization not available
    },
    'CASSCF': {
        'basis_function': '6-31G(d)',
        'memory_mb': 3000,
        'ncas': 4,                    # Number of active orbitals
        'nelecas': 4,                 # Number of active electrons
        'max_cycle_macro': 50,        # Orbital optimization iterations
        'max_cycle_micro': 3,         # CI solver iterations
        'natorb': True,               # Natural orbital transformation recommended
        'conv_tol': 1e-6,             # Energy convergence tolerance
        'conv_tol_grad': 1e-4,        # Gradient convergence tolerance
        'optimize_geometry': False    # Geometry optimization not available
    }
}

# Parameter validation constraints
# Defines min/max bounds, applicable methods, and UI constraints for each parameter
PARAMETER_CONSTRAINTS: Dict[str, Dict[str, Any]] = {
    'ncas': {
        'min': 1,
        'max': 20,
        'applicable_methods': ['CASCI', 'CASSCF'],
        'description': 'Number of active space orbitals (1-20)'
    },
    'nelecas': {
        'min': 1,
        'max': 40,
        'applicable_methods': ['CASCI', 'CASSCF'],
        'description': 'Number of active space electrons (1-40)'
    },
    'max_cycle_macro': {
        'min': 1,
        'max': 200,
        'applicable_methods': ['CASSCF'],
        'description': 'Maximum CASSCF macro iterations (1-200)'
    },
    'max_cycle_micro': {
        'min': 1,
        'max': 100,
        'applicable_methods': ['CASCI', 'CASSCF'],
        'description': 'Maximum CI micro iterations (1-100)'
    },
    'tddft_nstates': {
        'min': 1,
        'max': 50,
        'applicable_methods': ['TDDFT'],
        'description': 'Number of excited states to calculate (1-50)'
    },
    'optimize_geometry': {
        'disabled_for': ['TDDFT', 'CASCI', 'CASSCF', 'CCSD', 'CCSD_T'],
        'description': 'Geometry optimization is not available for these calculation methods'
    },
    'frozen_core': {
        'applicable_methods': ['CCSD', 'CCSD_T'],
        'description': 'Freeze core orbitals to reduce computational cost (recommended for CCSD/CCSD(T))'
    },
    'tddft_analyze_nto': {
        'applicable_methods': ['TDDFT'],
        'description': 'Perform Natural Transition Orbital analysis for excited states'
    },
    'cpu_cores': {
        'min': 1,
        'max': 32,
        'description': 'Number of CPU cores to use (1-32)'
    },
    'memory_mb': {
        'min': 128,
        'description': 'Memory allocation in megabytes (minimum 128 MB)'
    }
}


def get_method_defaults() -> Dict[str, Dict[str, Any]]:
    """Get all method default values.

    Returns:
        Dictionary mapping calculation method names to their default parameter values.

    Example:
        {
            'DFT': {'basis_function': '6-31G(d)', 'memory_mb': 2000, ...},
            'CCSD': {'basis_function': 'cc-pVDZ', 'memory_mb': 4000, ...},
            ...
        }
    """
    return METHOD_DEFAULTS


def get_parameter_constraints() -> Dict[str, Dict[str, Any]]:
    """Get all parameter constraints.

    Returns:
        Dictionary mapping parameter names to their constraints (min, max, applicable_methods, etc.).

    Example:
        {
            'ncas': {'min': 1, 'max': 20, 'applicable_methods': ['CASCI', 'CASSCF'], ...},
            'optimize_geometry': {'disabled_for': ['TDDFT', 'CASCI', ...], ...},
            ...
        }
    """
    return PARAMETER_CONSTRAINTS


def get_defaults_for_method(method: str) -> Dict[str, Any]:
    """Get default parameter values for a specific calculation method.

    Args:
        method: Calculation method name (e.g., 'DFT', 'CCSD', 'TDDFT')

    Returns:
        Dictionary of default parameter values for the method.
        Returns empty dict if method is not recognized.

    Example:
        >>> get_defaults_for_method('CCSD')
        {'basis_function': 'cc-pVDZ', 'memory_mb': 4000, 'frozen_core': True, ...}
    """
    return METHOD_DEFAULTS.get(method, {})


def is_parameter_applicable(param_name: str, method: str) -> bool:
    """Check if a parameter is applicable for a specific calculation method.

    Args:
        param_name: Name of the parameter (e.g., 'ncas', 'tddft_nstates')
        method: Calculation method name (e.g., 'CASCI', 'TDDFT')

    Returns:
        True if the parameter is applicable for the method, False otherwise.
        Returns True if no constraint exists (parameter is universally applicable).

    Example:
        >>> is_parameter_applicable('ncas', 'CASCI')
        True
        >>> is_parameter_applicable('ncas', 'DFT')
        False
    """
    constraint = PARAMETER_CONSTRAINTS.get(param_name)
    if not constraint:
        return True

    if 'applicable_methods' in constraint:
        return method in constraint['applicable_methods']

    return True


def is_parameter_disabled(param_name: str, method: str) -> bool:
    """Check if a parameter is disabled for a specific calculation method.

    Args:
        param_name: Name of the parameter (e.g., 'optimize_geometry')
        method: Calculation method name (e.g., 'TDDFT', 'CCSD')

    Returns:
        True if the parameter is disabled for the method, False otherwise.

    Example:
        >>> is_parameter_disabled('optimize_geometry', 'TDDFT')
        True
        >>> is_parameter_disabled('optimize_geometry', 'DFT')
        False
    """
    constraint = PARAMETER_CONSTRAINTS.get(param_name)
    if not constraint:
        return False

    if 'disabled_for' in constraint:
        return method in constraint['disabled_for']

    return False


def validate_parameter_value(param_name: str, value: Any) -> tuple[bool, str]:
    """Validate a parameter value against its constraints.

    Args:
        param_name: Name of the parameter
        value: Value to validate

    Returns:
        Tuple of (is_valid, error_message). error_message is empty string if valid.

    Example:
        >>> validate_parameter_value('ncas', 5)
        (True, '')
        >>> validate_parameter_value('ncas', 50)
        (False, 'Value 50 exceeds maximum of 20')
    """
    constraint = PARAMETER_CONSTRAINTS.get(param_name)
    if not constraint:
        return True, ''

    # Check minimum value
    if 'min' in constraint and value < constraint['min']:
        return False, f"Value {value} is below minimum of {constraint['min']}"

    # Check maximum value
    if 'max' in constraint and value > constraint['max']:
        return False, f"Value {value} exceeds maximum of {constraint['max']}"

    return True, ''
