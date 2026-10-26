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
        'optimize_geometry': True,
        'geomopt_maxsteps': 100,
        'geomopt_conv_energy': 1e-6,
        'density_fitting': False
    },
    'HF': {
        'basis_function': '6-31G(d)',
        'memory_mb': 2000,
        'optimize_geometry': True,
        'geomopt_maxsteps': 100,
        'geomopt_conv_energy': 1e-6,
        'density_fitting': False
    },
    'MP2': {
        'basis_function': '6-31G(d)',
        'memory_mb': 3000,
        'optimize_geometry': True,
        'geomopt_maxsteps': 100,
        'geomopt_conv_energy': 1e-6,
        'density_fitting': False
    },
    'CCSD': {
        'basis_function': 'cc-pVDZ',  # Correlation-consistent basis recommended for CCSD
        'memory_mb': 4000,            # Higher memory requirements
        'frozen_core': True,          # Frozen core approximation to reduce cost
        'density_fitting': False
    },
    'CCSD_T': {
        'basis_function': 'cc-pVDZ',  # Correlation-consistent basis recommended
        'memory_mb': 4000,            # Higher memory requirements
        'frozen_core': True,          # Frozen core approximation to reduce cost
        'density_fitting': False
    },
    'TDDFT': {
        'basis_function': '6-31G(d)',
        'exchange_correlation': 'B3LYP',
        'memory_mb': 2000,
        'tddft_nstates': 10,          # Default number of excited states
        'tddft_method': 'TDDFT',      # Full TDDFT (vs TDA approximation)
        'tddft_analyze_nto': False,   # NTO analysis off by default
        'density_fitting': False
    },
    'CASCI': {
        'basis_function': '6-31G(d)',
        'memory_mb': 3000,
        'ncas': 4,                    # Number of active orbitals
        'nelecas': 4,                 # Number of active electrons
        'natorb': True,               # Natural orbital transformation recommended
        'max_cycle_micro': 3          # CI solver iterations
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
        'conv_tol_grad': 1e-4         # Gradient convergence tolerance
    }
}

# Parameter validation constraints
# Defines min/max bounds, applicable methods, and UI constraints for each parameter
PARAMETER_CONSTRAINTS: Dict[str, Dict[str, Any]] = {
    'exchange_correlation': {
        'applicable_methods': ['DFT', 'TDDFT'],
        'description': 'Exchange-correlation functional (only applicable for DFT and TDDFT methods)'
    },
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
        'applicable_methods': ['DFT', 'HF', 'MP2'],
        'description': 'DFT, HF, and MP2 methods only'
    },
    'geomopt_maxsteps': {
        'min': 1,
        'max': 1000,
        'applicable_methods': ['DFT', 'HF', 'MP2'],
        'description': 'Maximum geometry optimization steps (1-1000)'
    },
    'geomopt_conv_energy': {
        'applicable_methods': ['DFT', 'HF', 'MP2'],
        'description': 'Energy convergence threshold in Hartree'
    },
    'frozen_core': {
        'applicable_methods': ['CCSD', 'CCSD_T'],
        'description': 'Freeze core orbitals to reduce computational cost (recommended for CCSD/CCSD(T))'
    },
    'tddft_analyze_nto': {
        'applicable_methods': ['TDDFT'],
        'description': 'Perform Natural Transition Orbital analysis for excited states'
    },
    'density_fitting': {
        'applicable_methods': ['DFT', 'HF', 'MP2', 'CCSD', 'CCSD_T', 'TDDFT'],
        'description': 'Enable density fitting (RI/DF) approximation for faster integral evaluation'
    },
    'auxiliary_basis': {
        'applicable_methods': ['DFT', 'HF', 'MP2', 'CCSD', 'CCSD_T', 'TDDFT'],
        'description': 'Auxiliary basis set for density fitting (auto-selected if not specified)'
    },
    'cpu_cores': {
        'min': 1,
        'max': 32,
        'description': 'Number of CPU cores to use (1-32)'
    },
    'memory_mb': {
        'min': 512,
        'max': 32768,
        'description': 'Memory allocation in megabytes (512-32768 MB)'
    },
    'charges': {
        'min': -10,
        'max': 10,
        'description': 'Molecular charge (-10 to 10)'
    },
    'spin': {
        'min': 0,
        'max': 10,
        'description': 'Spin (2S), number of unpaired electrons (0-10)'
    }
}

UNIVERSAL_PARAMS = {
    'xyz', 'calculation_method', 'basis_function', 'charges', 'spin',
    'solvent_method', 'solvent', 'name', 'cpu_cores', 'memory_mb',
    'ketcher_data', 'created_at'
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

    try:
        # Check minimum value
        if 'min' in constraint and value < constraint['min']:
            if constraint['min'] == 1:
                return False, f"Value {value} is below minimum of 1 and must be greater than 0"
            return False, f"Value {value} is below minimum of {constraint['min']}"

        # Check maximum value
        if 'max' in constraint and value > constraint['max']:
            return False, f"Value {value} exceeds maximum of {constraint['max']}"
    except TypeError:
        return False, f"Value {value} cannot be compared with numeric constraints"

    return True, ''


def validate_parameters_for_method(
    method: str,
    params: Dict[str, Any]
) -> tuple[bool, str]:
    """Validate that only applicable parameters are provided for a calculation method.

    This function performs strict validation by rejecting requests that contain
    parameters which are either:
    1. Not applicable to the specified method (e.g., 'ncas' for DFT)
    2. Explicitly disabled for the method (e.g., 'optimize_geometry' for TDDFT)

    Args:
        method: Calculation method name (e.g., 'DFT', 'CASCI', 'TDDFT')
        params: Dictionary of all provided parameters including their values

    Returns:
        Tuple of (is_valid, error_message). error_message is empty string if valid.

    Examples:
        >>> validate_parameters_for_method('DFT', {'xyz': '...', 'basis_function': '6-31G(d)'})
        (True, '')

        >>> validate_parameters_for_method('DFT', {'xyz': '...', 'ncas': 4})
        (False, "Parameter 'ncas' is not applicable for method 'DFT'. This parameter is only valid for: CASCI, CASSCF")

        >>> validate_parameters_for_method('TDDFT', {'xyz': '...', 'optimize_geometry': True})
        (False, "Parameter 'optimize_geometry' is disabled for method 'TDDFT'. Reason: Geometry optimization is not available for these calculation methods")
    """
    invalid_params = []
    invalid_values = []

    for param_name, param_value in params.items():
        # Skip None values (parameter not explicitly provided)
        if param_value is None:
            continue

        value_is_valid, value_error = validate_parameter_value(param_name, param_value)
        if not value_is_valid:
            invalid_values.append({
                'param': param_name,
                'error': value_error
            })

        # Skip universal parameters that apply to all methods
        if param_name in UNIVERSAL_PARAMS:
            continue

        # Check if parameter is applicable to this method
        if not is_parameter_applicable(param_name, method):
            constraint = PARAMETER_CONSTRAINTS.get(param_name, {})
            applicable_to = constraint.get('applicable_methods', [])
            invalid_params.append({
                'param': param_name,
                'value': param_value,
                'applicable_to': applicable_to
            })

    # Build comprehensive error message if any invalid parameters or values found
    if invalid_params or invalid_values:
        error_lines = []
        for invalid in invalid_values:
            error_lines.append(
                f"Parameter '{invalid['param']}' has invalid value: {invalid['error']}"
            )

        for invalid in invalid_params:
            param = invalid['param']
            applicable = ', '.join(invalid['applicable_to']) if invalid['applicable_to'] else 'none'
            error_lines.append(
                f"Parameter '{param}' is not applicable for method '{method}'. "
                f"This parameter is only valid for: {applicable}"
            )

        return False, '; '.join(error_lines)

    return True, ''
