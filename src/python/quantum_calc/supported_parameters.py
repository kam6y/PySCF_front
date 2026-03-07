"""Module for detecting supported quantum chemistry parameters from PySCF and other sources."""

import logging
from typing import Dict, List, Any
from pyscf import gto, dft
from data.solvent_properties import SOLVENT_DIELECTRIC
from quantum_calc.method_defaults import get_method_defaults, get_parameter_constraints

logger = logging.getLogger(__name__)


def get_supported_calculation_methods() -> List[str]:
    """Get list of supported quantum calculation methods."""
    return ['DFT', 'HF', 'MP2', 'CCSD', 'CCSD_T', 'TDDFT', 'CASCI', 'CASSCF']


def get_supported_basis_functions() -> Dict[str, List[str]]:
    """Get supported basis functions grouped by category."""
    return {
        'Minimal': [
            'STO-3G',
            'STO-6G',
            'MINAO'
        ],
        'Pople Style': [
            '3-21G',
            '4-31G',
            '6-31G',
            '6-31G(d)',
            '6-31G(d,p)',
            '6-31+G(d)',
            '6-31+G(d,p)',
            '6-31++G(d,p)',
            '6-311G',
            '6-311G(d)',
            '6-311G(d,p)',
            '6-311+G(d)',
            '6-311+G(d,p)',
            '6-311++G(d,p)'
        ],
        'Correlation Consistent': [
            'cc-pVDZ',
            'cc-pVTZ',
            'cc-pVQZ',
            'cc-pV5Z',
            'aug-cc-pVDZ',
            'aug-cc-pVTZ',
            'aug-cc-pVQZ',
            'aug-cc-pV5Z'
        ],
        'def2': [
            'def2-SVP',
            'def2-SVPD',
            'def2-TZVP',
            'def2-TZVPP',
            'def2-TZVPD',
            'def2-QZVP',
            'def2-QZVPP',
            'def2-QZVPD'
        ]
    }


def get_supported_auxiliary_basis_functions() -> Dict[str, List[str]]:
    """Get supported auxiliary basis functions for density fitting, grouped by category."""
    return {
        'def2 RI': [
            'def2-svp-ri',
            'def2-tzvp-ri',
            'def2-tzvpp-ri',
            'def2-qzvp-ri',
            'def2-qzvpp-ri',
        ],
        'cc RI': [
            'cc-pvdz-ri',
            'cc-pvtz-ri',
            'cc-pvqz-ri',
            'aug-cc-pvdz-ri',
            'aug-cc-pvtz-ri',
            'aug-cc-pvqz-ri',
        ],
        'Weigend': [
            'weigend',
        ],
    }


def get_supported_exchange_correlation() -> Dict[str, List[str]]:
    """Get supported exchange-correlation functionals grouped by category."""
    return {
        'LDA': [
            'SVWN'
        ],
        'GGA': [
            'PBE',
            'BLYP',
            'BP86',
            'PW91',
            'revPBE',
            'PBEsol',
            'RPBE',
            'OLYP',
            'BPBE'
        ],
        'Meta-GGA': [
            'M06-L',
            'TPSS',
            'revTPSS',
            'SCAN',
            'rSCAN',
            'r2SCAN'
        ],
        'Hybrid (GGA)': [
            'B3LYP',
            'B3PW91',
            'PBE0',
            'B3P86',
            'O3LYP'
        ],
        'Hybrid (Meta-GGA)': [
            'M06',
            'M06-2X',
            'TPSSh'
        ],
        'Range-separated Hybrid': [
            'CAM-B3LYP',
            'wB97XD'
        ]
    }


def get_supported_solvent_methods() -> List[str]:
    """Get list of supported solvent effect methods."""
    return ['none', 'ief-pcm', 'c-pcm', 'cosmo', 'ssvpe', 'ddcosmo']


def get_supported_solvents() -> Dict[str, List[Dict[str, Any]]]:
    """Get supported solvents grouped by category with their properties."""
    return {
        'Highly Polar': [
            {'value': 'water', 'display': 'Water (78.36)', 'dielectric_constant': 78.3553},
            {'value': 'dimethylsulfoxide', 'display': 'Dimethylsulfoxide (46.83)', 'dielectric_constant': 46.826},
            {'value': 'n,n-dimethylformamide', 'display': 'N,N-Dimethylformamide (37.22)', 'dielectric_constant': 37.219},
            {'value': 'nitromethane', 'display': 'Nitromethane (36.56)', 'dielectric_constant': 36.562}
        ],
        'Protic Solvents': [
            {'value': 'methanol', 'display': 'Methanol (32.61)', 'dielectric_constant': 32.613},
            {'value': 'ethanol', 'display': 'Ethanol (24.85)', 'dielectric_constant': 24.852}
        ],
        'Polar Aprotic': [
            {'value': 'acetone', 'display': 'Acetone (20.49)', 'dielectric_constant': 20.493},
            {'value': 'dichloroethane', 'display': 'Dichloroethane (10.13)', 'dielectric_constant': 10.125},
            {'value': 'dichloromethane', 'display': 'Dichloromethane (8.93)', 'dielectric_constant': 8.93},
            {'value': 'tetrahydrofuran', 'display': 'Tetrahydrofuran (7.43)', 'dielectric_constant': 7.4297},
            {'value': 'chlorobenzene', 'display': 'Chlorobenzene (5.70)', 'dielectric_constant': 5.6968}
        ],
        'Moderately Polar': [
            {'value': 'chloroform', 'display': 'Chloroform (4.71)', 'dielectric_constant': 4.7113},
            {'value': 'diethylether', 'display': 'Diethylether (4.24)', 'dielectric_constant': 4.2400}
        ],
        'Nonpolar': [
            {'value': 'toluene', 'display': 'Toluene (2.37)', 'dielectric_constant': 2.3741},
            {'value': 'benzene', 'display': 'Benzene (2.27)', 'dielectric_constant': 2.2706},
            {'value': '1,4-dioxane', 'display': '1,4-Dioxane (2.21)', 'dielectric_constant': 2.2099},
            {'value': 'cyclohexane', 'display': 'Cyclohexane (2.02)', 'dielectric_constant': 2.0160}
        ]
    }


def get_supported_tddft_methods() -> List[str]:
    """Get list of supported TDDFT calculation methods."""
    return ['TDDFT', 'TDA']


def validate_basis_function(basis: str) -> bool:
    """Validate if a basis function is supported by PySCF."""
    try:
        mol = gto.Mole()
        mol.atom = 'H 0 0 0; H 0 0 0.74'
        mol.basis = basis
        mol.build(verbose=0)
        return True
    except Exception as e:
        logger.debug(f"Basis function {basis} not supported: {e}")
        return False


def validate_exchange_correlation(xc: str) -> bool:
    """Validate if an exchange-correlation functional is supported by PySCF."""
    try:
        mol = gto.Mole()
        mol.atom = 'H 0 0 0; H 0 0 0.74'
        mol.basis = 'sto-3g'
        mol.build(verbose=0)

        mf = dft.RKS(mol)
        mf.xc = xc
        mf.max_cycle = 1
        mf.verbose = 0
        mf.kernel()
        return True
    except Exception as e:
        logger.debug(f"Exchange-correlation functional {xc} not supported: {e}")
        return False


def get_all_supported_parameters() -> Dict[str, Any]:
    """Get all supported parameters in a structured format."""
    return {
        'calculation_methods': get_supported_calculation_methods(),
        'basis_functions': get_supported_basis_functions(),
        'auxiliary_basis_functions': get_supported_auxiliary_basis_functions(),
        'exchange_correlation': get_supported_exchange_correlation(),
        'solvent_methods': get_supported_solvent_methods(),
        'solvents': get_supported_solvents(),
        'tddft_methods': get_supported_tddft_methods(),
        'method_defaults': get_method_defaults(),
        'parameter_constraints': get_parameter_constraints()
    }
