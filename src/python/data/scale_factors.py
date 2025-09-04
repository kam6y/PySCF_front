"""Scale factors for vibrational frequency analysis in quantum chemistry calculations.

This module provides scale factors for correcting computed vibrational frequencies
to match experimental values more closely. The scale factors are method- and 
basis set-specific, based on empirical comparisons with experimental data.

References:
- NIST Computational Chemistry Comparison and Benchmark Database
- https://cccbdb.nist.gov/vibscalejust.asp
"""

from typing import Dict, Optional, Tuple
import logging

logger = logging.getLogger(__name__)

# Scale factors for different computational methods and basis sets
# Format: method -> basis_set -> scale_factor
SCALE_FACTORS: Dict[str, Dict[str, float]] = {
    'HF': {
        '6-31G*': 0.899,
        '6-31G**': 0.903,
        '6-31+G**': 0.904,
        'cc-pVDZ': 0.908,
        'cc-pVTZ': 0.910
    },
    'B3LYP': {
        '6-31G*': 0.960,
        '6-31G**': 0.961,
        '6-31+G**': 0.964,
        'cc-pVDZ': 0.970,
        'cc-pVTZ': 0.967
    },
    'M06-2X': {
        '6-31G*': 0.947,
        '6-31G**': 0.950,
        '6-31+G**': 0.952,
        'cc-pVTZ': 0.955
    },
    'PBEPBE': {
        '6-31G*': 0.986,
        '6-31G**': 0.986,
        '6-31+G**': 0.989,
        'cc-pVDZ': 0.994,
        'cc-pVTZ': 0.993
    },
    'TPSSh': {
        '6-31G*': 0.959,
        '6-31G**': 0.959,
        '6-31+G**': 0.963,
        'cc-pVDZ': 0.972,
        'cc-pVTZ': 0.968
    },
    'wB97X-D': {
        '6-31G*': 0.949,
        '6-31+G**': 0.952,
        'cc-pVDZ': 0.953,
        'cc-pVTZ': 0.956
    },
    'MP2': {
        '6-31G*': 0.943,
        '6-31G**': 0.937,
        '6-31+G**': 0.941,
        'cc-pVDZ': 0.953,
        'cc-pVTZ': 0.950
    },
    'CID': {  # Configuration Interaction with Double excitations
        '6-31G*': 0.924,
        '6-31G**': 0.924,
        '6-31+G**': 0.924,
        'cc-pVDZ': 0.924,
        'cc-pVTZ': 0.927
    }
}

# Alternative method name mappings for PySCF compatibility
METHOD_ALIASES: Dict[str, str] = {
    'B3LYP': 'B3LYP',
    'b3lyp': 'B3LYP',
    'HF': 'HF',
    'hf': 'HF',
    'RHF': 'HF',
    'rhf': 'HF',
    'UHF': 'HF',
    'uhf': 'HF',
    'MP2': 'MP2',
    'mp2': 'MP2',
    'PBE': 'PBEPBE',
    'pbe': 'PBEPBE',
    'PBEPBE': 'PBEPBE',
    'pbepbe': 'PBEPBE',
    'M062X': 'M06-2X',
    'm062x': 'M06-2X',
    'M06-2X': 'M06-2X',
    'm06-2x': 'M06-2X',
    'TPSSh': 'TPSSh',
    'tpssh': 'TPSSh',
    'wB97XD': 'wB97X-D',
    'wb97xd': 'wB97X-D',
    'wB97X-D': 'wB97X-D',
    'wb97x-d': 'wB97X-D',
    'CCSD': 'CID',  # Approximate mapping for CCSD to CID
    'ccsd': 'CID'
}

# Basis set name mappings for consistency
BASIS_ALIASES: Dict[str, str] = {
    '6-31G*': '6-31G*',
    '6-31G(d)': '6-31G*',
    '6-31g*': '6-31G*',
    '6-31g(d)': '6-31G*',
    '6-31G**': '6-31G**',
    '6-31G(d,p)': '6-31G**',
    '6-31g**': '6-31G**',
    '6-31g(d,p)': '6-31G**',
    '6-31+G**': '6-31+G**',
    '6-31+G(d,p)': '6-31+G**',
    '6-31+g**': '6-31+G**',
    '6-31+g(d,p)': '6-31+G**',
    'cc-pVDZ': 'cc-pVDZ',
    'cc-pvdz': 'cc-pVDZ',
    'cc-pVTZ': 'cc-pVTZ',
    'cc-pvtz': 'cc-pVTZ'
}


def get_scale_factor(method: str, basis_set: str) -> Optional[float]:
    """
    Get the scale factor for a given computational method and basis set.
    
    Args:
        method: Computational method (e.g., 'B3LYP', 'HF', 'MP2')
        basis_set: Basis set (e.g., '6-31G*', 'cc-pVDZ')
        
    Returns:
        Scale factor as float, or None if not available
        
    Examples:
        >>> get_scale_factor('B3LYP', '6-31G*')
        0.960
        >>> get_scale_factor('HF', 'cc-pVDZ') 
        0.908
        >>> get_scale_factor('unknown', '6-31G*')
        None
    """
    # Normalize method name
    normalized_method = METHOD_ALIASES.get(method, method)
    
    # Normalize basis set name
    normalized_basis = BASIS_ALIASES.get(basis_set, basis_set)
    
    # Look up scale factor
    if normalized_method in SCALE_FACTORS:
        method_data = SCALE_FACTORS[normalized_method]
        scale_factor = method_data.get(normalized_basis)
        
        if scale_factor is not None:
            logger.debug(f"Found scale factor {scale_factor} for {normalized_method}/{normalized_basis}")
            return scale_factor
        else:
            logger.warning(f"No scale factor available for basis set '{normalized_basis}' with method '{normalized_method}'")
    else:
        logger.warning(f"No scale factors available for method '{normalized_method}'")
    
    return None


def get_recommended_scale_factor(method: str, basis_set: str) -> Tuple[float, str]:
    """
    Get a scale factor for the given method/basis, with fallback recommendations.
    
    Args:
        method: Computational method
        basis_set: Basis set
        
    Returns:
        Tuple of (scale_factor, recommendation_message)
        
    Examples:
        >>> get_recommended_scale_factor('B3LYP', '6-31G*')
        (0.960, 'Exact match found')
        >>> get_recommended_scale_factor('B3LYP', 'STO-3G')
        (0.960, 'Using fallback: 6-31G* scale factor for B3LYP')
    """
    # Try exact match first
    scale_factor = get_scale_factor(method, basis_set)
    if scale_factor is not None:
        return scale_factor, "Exact match found"
    
    # Normalize method name for fallback logic
    normalized_method = METHOD_ALIASES.get(method, method)
    
    # Fallback strategy 1: Use the most common basis set for this method
    if normalized_method in SCALE_FACTORS:
        method_data = SCALE_FACTORS[normalized_method]
        
        # Preference order for fallback basis sets
        fallback_order = ['6-31G*', '6-31G**', 'cc-pVDZ', '6-31+G**', 'cc-pVTZ']
        
        for fallback_basis in fallback_order:
            if fallback_basis in method_data:
                fallback_factor = method_data[fallback_basis]
                message = f"Using fallback: {fallback_basis} scale factor for {normalized_method}"
                logger.info(message)
                return fallback_factor, message
    
    # Fallback strategy 2: Use B3LYP scale factor as general approximation
    if normalized_method != 'B3LYP':
        b3lyp_factor = get_scale_factor('B3LYP', basis_set)
        if b3lyp_factor is not None:
            message = f"Using B3LYP scale factor as approximation for {normalized_method}"
            logger.info(message)
            return b3lyp_factor, message
        
        # Use B3LYP with 6-31G* as ultimate fallback
        b3lyp_631g_factor = get_scale_factor('B3LYP', '6-31G*')
        if b3lyp_631g_factor is not None:
            message = f"Using B3LYP/6-31G* scale factor ({b3lyp_631g_factor}) as general approximation"
            logger.warning(message)
            return b3lyp_631g_factor, message
    
    # Ultimate fallback: use 0.96 (typical for DFT methods)
    default_factor = 0.96
    message = f"No specific scale factor found. Using default value ({default_factor})"
    logger.warning(message)
    return default_factor, message


def list_available_methods() -> Dict[str, list]:
    """
    List all available methods and their supported basis sets.
    
    Returns:
        Dictionary mapping method names to lists of supported basis sets
    """
    return {method: list(basis_data.keys()) for method, basis_data in SCALE_FACTORS.items()}


def validate_method_basis_combination(method: str, basis_set: str) -> Tuple[bool, str]:
    """
    Validate if a method/basis combination has a scale factor available.
    
    Args:
        method: Computational method
        basis_set: Basis set
        
    Returns:
        Tuple of (is_available, message)
    """
    scale_factor = get_scale_factor(method, basis_set)
    
    if scale_factor is not None:
        return True, f"Scale factor available: {scale_factor}"
    else:
        # Check if method exists with other basis sets
        normalized_method = METHOD_ALIASES.get(method, method)
        if normalized_method in SCALE_FACTORS:
            available_basis = list(SCALE_FACTORS[normalized_method].keys())
            return False, f"Method '{normalized_method}' is supported, but not with basis '{basis_set}'. Available basis sets: {available_basis}"
        else:
            available_methods = list(SCALE_FACTORS.keys())
            return False, f"Method '{method}' not supported. Available methods: {available_methods}"