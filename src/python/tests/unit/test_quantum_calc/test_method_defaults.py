"""Unit tests for method defaults and parameter constraints."""

import pytest
from quantum_calc.method_defaults import (
    get_method_defaults,
    get_parameter_constraints,
    get_defaults_for_method,
    is_parameter_applicable,
    is_parameter_disabled,
    validate_parameter_value,
    METHOD_DEFAULTS,
    PARAMETER_CONSTRAINTS,
)


class TestMethodDefaults:
    """Tests for method default values."""

    def test_get_method_defaults_returns_dict(self):
        """Test that get_method_defaults returns a dictionary."""
        defaults = get_method_defaults()
        assert isinstance(defaults, dict)
        assert len(defaults) > 0

    def test_all_calculation_methods_have_defaults(self):
        """Test that all supported calculation methods have defaults."""
        expected_methods = ['DFT', 'HF', 'MP2', 'CCSD', 'CCSD_T', 'TDDFT', 'CASCI', 'CASSCF']
        defaults = get_method_defaults()

        for method in expected_methods:
            assert method in defaults, f"Method {method} missing defaults"
            assert isinstance(defaults[method], dict), f"Defaults for {method} should be a dict"

    def test_dft_defaults(self):
        """Test DFT default values."""
        defaults = get_defaults_for_method('DFT')
        assert defaults['basis_function'] == '6-31G(d)'
        assert defaults['exchange_correlation'] == 'B3LYP'
        assert defaults['memory_mb'] == 2000
        assert defaults['optimize_geometry'] is True

    def test_ccsd_defaults(self):
        """Test CCSD default values."""
        defaults = get_defaults_for_method('CCSD')
        assert defaults['basis_function'] == 'cc-pVDZ'
        assert defaults['memory_mb'] == 4000
        assert defaults['frozen_core'] is True
        assert defaults['optimize_geometry'] is False

    def test_ccsd_t_defaults(self):
        """Test CCSD(T) default values."""
        defaults = get_defaults_for_method('CCSD_T')
        assert defaults['basis_function'] == 'cc-pVDZ'
        assert defaults['memory_mb'] == 4000
        assert defaults['frozen_core'] is True
        assert defaults['optimize_geometry'] is False

    def test_tddft_defaults(self):
        """Test TDDFT default values."""
        defaults = get_defaults_for_method('TDDFT')
        assert defaults['basis_function'] == '6-31G(d)'
        assert defaults['exchange_correlation'] == 'B3LYP'
        assert defaults['memory_mb'] == 2000
        assert defaults['tddft_nstates'] == 10
        assert defaults['tddft_method'] == 'TDDFT'
        assert defaults['tddft_analyze_nto'] is False
        assert defaults['optimize_geometry'] is False

    def test_casci_defaults(self):
        """Test CASCI default values."""
        defaults = get_defaults_for_method('CASCI')
        assert defaults['basis_function'] == '6-31G(d)'
        assert defaults['memory_mb'] == 3000
        assert defaults['ncas'] == 4
        assert defaults['nelecas'] == 4
        assert defaults['natorb'] is True
        assert defaults['max_cycle_micro'] == 3
        assert defaults['optimize_geometry'] is False

    def test_casscf_defaults(self):
        """Test CASSCF default values."""
        defaults = get_defaults_for_method('CASSCF')
        assert defaults['basis_function'] == '6-31G(d)'
        assert defaults['memory_mb'] == 3000
        assert defaults['ncas'] == 4
        assert defaults['nelecas'] == 4
        assert defaults['max_cycle_macro'] == 50
        assert defaults['max_cycle_micro'] == 3
        assert defaults['natorb'] is True
        assert defaults['conv_tol'] == 1e-6
        assert defaults['conv_tol_grad'] == 1e-4
        assert defaults['optimize_geometry'] is False

    def test_hf_defaults(self):
        """Test HF default values."""
        defaults = get_defaults_for_method('HF')
        assert defaults['basis_function'] == '6-31G(d)'
        assert defaults['memory_mb'] == 2000
        assert defaults['optimize_geometry'] is True

    def test_mp2_defaults(self):
        """Test MP2 default values."""
        defaults = get_defaults_for_method('MP2')
        assert defaults['basis_function'] == '6-31G(d)'
        assert defaults['memory_mb'] == 3000
        assert defaults['optimize_geometry'] is True

    def test_get_defaults_for_unknown_method(self):
        """Test getting defaults for an unknown method returns empty dict."""
        defaults = get_defaults_for_method('UNKNOWN_METHOD')
        assert defaults == {}


class TestParameterConstraints:
    """Tests for parameter constraints."""

    def test_get_parameter_constraints_returns_dict(self):
        """Test that get_parameter_constraints returns a dictionary."""
        constraints = get_parameter_constraints()
        assert isinstance(constraints, dict)
        assert len(constraints) > 0

    def test_ncas_constraint(self):
        """Test ncas parameter constraint."""
        constraints = get_parameter_constraints()
        assert 'ncas' in constraints
        assert constraints['ncas']['min'] == 1
        assert constraints['ncas']['max'] == 20
        assert 'CASCI' in constraints['ncas']['applicable_methods']
        assert 'CASSCF' in constraints['ncas']['applicable_methods']

    def test_nelecas_constraint(self):
        """Test nelecas parameter constraint."""
        constraints = get_parameter_constraints()
        assert 'nelecas' in constraints
        assert constraints['nelecas']['min'] == 1
        assert constraints['nelecas']['max'] == 40
        assert 'CASCI' in constraints['nelecas']['applicable_methods']
        assert 'CASSCF' in constraints['nelecas']['applicable_methods']

    def test_optimize_geometry_constraint(self):
        """Test optimize_geometry parameter constraint."""
        constraints = get_parameter_constraints()
        assert 'optimize_geometry' in constraints
        disabled_methods = constraints['optimize_geometry']['disabled_for']
        assert 'TDDFT' in disabled_methods
        assert 'CASCI' in disabled_methods
        assert 'CASSCF' in disabled_methods
        assert 'CCSD' in disabled_methods
        assert 'CCSD_T' in disabled_methods

    def test_tddft_nstates_constraint(self):
        """Test tddft_nstates parameter constraint."""
        constraints = get_parameter_constraints()
        assert 'tddft_nstates' in constraints
        assert constraints['tddft_nstates']['min'] == 1
        assert constraints['tddft_nstates']['max'] == 50
        assert 'TDDFT' in constraints['tddft_nstates']['applicable_methods']

    def test_cpu_cores_constraint(self):
        """Test cpu_cores parameter constraint."""
        constraints = get_parameter_constraints()
        assert 'cpu_cores' in constraints
        assert constraints['cpu_cores']['min'] == 1
        assert constraints['cpu_cores']['max'] == 32

    def test_memory_mb_constraint(self):
        """Test memory_mb parameter constraint."""
        constraints = get_parameter_constraints()
        assert 'memory_mb' in constraints
        assert constraints['memory_mb']['min'] == 128


class TestParameterApplicability:
    """Tests for parameter applicability checks."""

    def test_ncas_applicable_to_casci(self):
        """Test that ncas is applicable to CASCI."""
        assert is_parameter_applicable('ncas', 'CASCI') is True

    def test_ncas_not_applicable_to_dft(self):
        """Test that ncas is not applicable to DFT."""
        assert is_parameter_applicable('ncas', 'DFT') is False

    def test_tddft_nstates_applicable_to_tddft(self):
        """Test that tddft_nstates is applicable to TDDFT."""
        assert is_parameter_applicable('tddft_nstates', 'TDDFT') is True

    def test_tddft_nstates_not_applicable_to_dft(self):
        """Test that tddft_nstates is not applicable to DFT."""
        assert is_parameter_applicable('tddft_nstates', 'DFT') is False

    def test_universal_parameter_is_always_applicable(self):
        """Test that parameters without applicability constraints are always applicable."""
        # cpu_cores and memory_mb have no applicability constraints
        assert is_parameter_applicable('cpu_cores', 'DFT') is True
        assert is_parameter_applicable('cpu_cores', 'CCSD') is True
        assert is_parameter_applicable('memory_mb', 'TDDFT') is True


class TestParameterDisabled:
    """Tests for parameter disabled checks."""

    def test_optimize_geometry_disabled_for_tddft(self):
        """Test that optimize_geometry is disabled for TDDFT."""
        assert is_parameter_disabled('optimize_geometry', 'TDDFT') is True

    def test_optimize_geometry_disabled_for_casci(self):
        """Test that optimize_geometry is disabled for CASCI."""
        assert is_parameter_disabled('optimize_geometry', 'CASCI') is True

    def test_optimize_geometry_not_disabled_for_dft(self):
        """Test that optimize_geometry is not disabled for DFT."""
        assert is_parameter_disabled('optimize_geometry', 'DFT') is False

    def test_optimize_geometry_disabled_for_ccsd(self):
        """Test that optimize_geometry is disabled for CCSD."""
        assert is_parameter_disabled('optimize_geometry', 'CCSD') is True

    def test_parameter_without_disabled_constraint(self):
        """Test that parameters without disabled constraints are never disabled."""
        assert is_parameter_disabled('ncas', 'CASCI') is False
        assert is_parameter_disabled('cpu_cores', 'DFT') is False


class TestParameterValidation:
    """Tests for parameter value validation."""

    def test_validate_ncas_within_bounds(self):
        """Test validating ncas within valid bounds."""
        is_valid, error = validate_parameter_value('ncas', 10)
        assert is_valid is True
        assert error == ''

    def test_validate_ncas_below_min(self):
        """Test validating ncas below minimum."""
        is_valid, error = validate_parameter_value('ncas', 0)
        assert is_valid is False
        assert 'below minimum' in error

    def test_validate_ncas_above_max(self):
        """Test validating ncas above maximum."""
        is_valid, error = validate_parameter_value('ncas', 25)
        assert is_valid is False
        assert 'exceeds maximum' in error

    def test_validate_nelecas_within_bounds(self):
        """Test validating nelecas within valid bounds."""
        is_valid, error = validate_parameter_value('nelecas', 20)
        assert is_valid is True
        assert error == ''

    def test_validate_cpu_cores_within_bounds(self):
        """Test validating cpu_cores within valid bounds."""
        is_valid, error = validate_parameter_value('cpu_cores', 8)
        assert is_valid is True
        assert error == ''

    def test_validate_cpu_cores_above_max(self):
        """Test validating cpu_cores above maximum."""
        is_valid, error = validate_parameter_value('cpu_cores', 64)
        assert is_valid is False
        assert 'exceeds maximum' in error

    def test_validate_memory_mb_within_bounds(self):
        """Test validating memory_mb within valid bounds."""
        is_valid, error = validate_parameter_value('memory_mb', 2000)
        assert is_valid is True
        assert error == ''

    def test_validate_memory_mb_below_min(self):
        """Test validating memory_mb below minimum."""
        is_valid, error = validate_parameter_value('memory_mb', 64)
        assert is_valid is False
        assert 'below minimum' in error

    def test_validate_unknown_parameter(self):
        """Test validating an unknown parameter (should always pass)."""
        is_valid, error = validate_parameter_value('unknown_param', 999)
        assert is_valid is True
        assert error == ''


class TestDataIntegrity:
    """Tests for data integrity and consistency."""

    def test_all_methods_in_defaults_dict(self):
        """Test that METHOD_DEFAULTS contains all expected methods."""
        expected_methods = {'DFT', 'HF', 'MP2', 'CCSD', 'CCSD_T', 'TDDFT', 'CASCI', 'CASSCF'}
        actual_methods = set(METHOD_DEFAULTS.keys())
        assert actual_methods == expected_methods

    def test_all_constraints_have_description(self):
        """Test that all constraints have a description field."""
        for param_name, constraint in PARAMETER_CONSTRAINTS.items():
            assert 'description' in constraint, f"Parameter {param_name} missing description"
            assert isinstance(constraint['description'], str)
            assert len(constraint['description']) > 0

    def test_method_defaults_basis_function(self):
        """Test that all methods specify a basis_function."""
        for method, defaults in METHOD_DEFAULTS.items():
            assert 'basis_function' in defaults, f"Method {method} missing basis_function"
            assert isinstance(defaults['basis_function'], str)

    def test_method_defaults_memory_mb(self):
        """Test that all methods specify memory_mb."""
        for method, defaults in METHOD_DEFAULTS.items():
            assert 'memory_mb' in defaults, f"Method {method} missing memory_mb"
            assert isinstance(defaults['memory_mb'], int)
            assert defaults['memory_mb'] >= 128

    def test_method_defaults_optimize_geometry(self):
        """Test that all methods specify optimize_geometry."""
        for method, defaults in METHOD_DEFAULTS.items():
            assert 'optimize_geometry' in defaults, f"Method {method} missing optimize_geometry"
            assert isinstance(defaults['optimize_geometry'], bool)

    def test_dft_tddft_have_exchange_correlation(self):
        """Test that DFT and TDDFT specify exchange_correlation."""
        assert 'exchange_correlation' in METHOD_DEFAULTS['DFT']
        assert 'exchange_correlation' in METHOD_DEFAULTS['TDDFT']
        assert METHOD_DEFAULTS['DFT']['exchange_correlation'] == 'B3LYP'
        assert METHOD_DEFAULTS['TDDFT']['exchange_correlation'] == 'B3LYP'

    def test_constraint_min_less_than_max(self):
        """Test that min is always less than max in constraints."""
        for param_name, constraint in PARAMETER_CONSTRAINTS.items():
            if 'min' in constraint and 'max' in constraint:
                assert constraint['min'] < constraint['max'], \
                    f"Parameter {param_name}: min must be less than max"
