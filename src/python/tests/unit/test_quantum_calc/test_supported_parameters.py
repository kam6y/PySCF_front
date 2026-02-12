"""
Validation tests for supported PySCF parameters.

Tests that all basis sets and XC functionals listed in supported_parameters.py
actually work with PySCF. Any parameters that fail should be removed from
the supported lists.
"""

import pytest
from pyscf import gto, dft
from quantum_calc.supported_parameters import (
    get_supported_basis_functions,
    get_supported_exchange_correlation,
)

# Flatten all basis sets from categories
ALL_BASIS_SETS = [
    basis
    for bases in get_supported_basis_functions().values()
    for basis in bases
]

# Flatten all XC functionals from categories
ALL_XC_FUNCTIONALS = [
    xc
    for xcs in get_supported_exchange_correlation().values()
    for xc in xcs
]

H2_ATOM = 'H 0 0 0; H 0 0 0.74'


@pytest.mark.parametrize("basis", ALL_BASIS_SETS)
def test_basis_function_is_valid(basis):
    """
    GIVEN a basis set from supported_parameters.py
    WHEN we build an H2 molecule with that basis
    THEN PySCF should successfully build the molecule
    """
    mol = gto.Mole()
    mol.atom = H2_ATOM
    mol.basis = basis
    mol.build(verbose=0)


@pytest.mark.parametrize("xc", ALL_XC_FUNCTIONALS)
def test_xc_functional_is_valid(xc):
    """
    GIVEN an XC functional from supported_parameters.py
    WHEN we run one SCF iteration on H2/STO-3G with that functional
    THEN PySCF should execute without errors
    """
    mol = gto.Mole()
    mol.atom = H2_ATOM
    mol.basis = 'sto-3g'
    mol.build(verbose=0)

    mf = dft.RKS(mol)
    mf.xc = xc
    mf.max_cycle = 1
    mf.verbose = 0
    mf.kernel()
