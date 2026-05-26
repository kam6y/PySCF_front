"""Unit tests for method defaults and parameter constraints."""

import pytest

from quantum_calc.method_defaults import (
    get_method_defaults,
    get_parameter_constraints,
    get_defaults_for_method,
    is_parameter_applicable,
    validate_parameter_value,
    validate_parameters_for_method,
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

    @pytest.mark.parametrize(
        ("method", "expected_values", "absent_keys"),
        [
            (
                "DFT",
                {
                    "basis_function": "6-31G(d)",
                    "exchange_correlation": "B3LYP",
                    "memory_mb": 2000,
                    "optimize_geometry": True,
                },
                (),
            ),
            (
                "HF",
                {
                    "basis_function": "6-31G(d)",
                    "memory_mb": 2000,
                    "optimize_geometry": True,
                },
                (),
            ),
            (
                "MP2",
                {
                    "basis_function": "6-31G(d)",
                    "memory_mb": 3000,
                    "optimize_geometry": True,
                },
                (),
            ),
            (
                "CCSD",
                {
                    "basis_function": "cc-pVDZ",
                    "memory_mb": 4000,
                    "frozen_core": True,
                },
                ("optimize_geometry",),
            ),
            (
                "CCSD_T",
                {
                    "basis_function": "cc-pVDZ",
                    "memory_mb": 4000,
                    "frozen_core": True,
                },
                ("optimize_geometry",),
            ),
            (
                "TDDFT",
                {
                    "basis_function": "6-31G(d)",
                    "exchange_correlation": "B3LYP",
                    "memory_mb": 2000,
                    "tddft_nstates": 10,
                    "tddft_method": "TDDFT",
                    "tddft_analyze_nto": False,
                },
                ("optimize_geometry",),
            ),
            (
                "CASCI",
                {
                    "basis_function": "6-31G(d)",
                    "memory_mb": 3000,
                    "ncas": 4,
                    "nelecas": 4,
                    "natorb": True,
                    "max_cycle_micro": 3,
                },
                ("optimize_geometry",),
            ),
            (
                "CASSCF",
                {
                    "basis_function": "6-31G(d)",
                    "memory_mb": 3000,
                    "ncas": 4,
                    "nelecas": 4,
                    "max_cycle_macro": 50,
                    "max_cycle_micro": 3,
                    "natorb": True,
                    "conv_tol": 1e-6,
                    "conv_tol_grad": 1e-4,
                },
                ("optimize_geometry",),
            ),
        ],
    )
    def test_method_specific_defaults(
        self,
        method,
        expected_values,
        absent_keys,
    ):
        """Test representative method-specific default values."""
        defaults = get_defaults_for_method(method)

        for key, expected_value in expected_values.items():
            assert defaults[key] == expected_value

        for key in absent_keys:
            assert key not in defaults

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

    @pytest.mark.parametrize(
        ("param", "expected_values"),
        [
            ("ncas", {"min": 1, "max": 20, "applicable_methods": ["CASCI", "CASSCF"]}),
            ("nelecas", {"min": 1, "max": 40, "applicable_methods": ["CASCI", "CASSCF"]}),
            ("optimize_geometry", {"applicable_methods": ["DFT", "HF", "MP2"]}),
            ("tddft_nstates", {"min": 1, "max": 50, "applicable_methods": ["TDDFT"]}),
            ("cpu_cores", {"min": 1, "max": 32}),
            ("memory_mb", {"min": 512, "max": 32768}),
            ("charges", {"min": -10, "max": 10}),
            ("spin", {"min": 0, "max": 10}),
        ],
    )
    def test_parameter_constraints(self, param, expected_values):
        """Test representative parameter constraints."""
        constraints = get_parameter_constraints()

        assert param in constraints
        for key, expected_value in expected_values.items():
            assert constraints[param][key] == expected_value


class TestParameterApplicability:
    """Tests for parameter applicability checks."""

    @pytest.mark.parametrize(
        ("param", "method", "expected"),
        [
            ("ncas", "CASCI", True),
            ("ncas", "DFT", False),
            ("tddft_nstates", "TDDFT", True),
            ("tddft_nstates", "DFT", False),
            ("exchange_correlation", "DFT", True),
            ("exchange_correlation", "TDDFT", True),
            ("exchange_correlation", "HF", False),
            ("exchange_correlation", "MP2", False),
            ("cpu_cores", "DFT", True),
            ("cpu_cores", "CCSD", True),
            ("memory_mb", "TDDFT", True),
            ("optimize_geometry", "DFT", True),
            ("optimize_geometry", "HF", True),
            ("optimize_geometry", "MP2", True),
            ("optimize_geometry", "TDDFT", False),
            ("optimize_geometry", "CASCI", False),
            ("optimize_geometry", "CCSD", False),
        ],
    )
    def test_parameter_applicability(self, param, method, expected):
        """Test method-specific and universal parameter applicability."""
        assert is_parameter_applicable(param, method) is expected


class TestParameterValidation:
    """Tests for parameter value validation."""

    @pytest.mark.parametrize(
        ("param", "value", "expected_valid", "expected_error"),
        [
            ("ncas", 10, True, ""),
            ("ncas", 0, False, "below minimum"),
            ("ncas", 25, False, "exceeds maximum"),
            ("nelecas", 20, True, ""),
            ("cpu_cores", 8, True, ""),
            ("cpu_cores", 1, True, ""),
            ("cpu_cores", 32, True, ""),
            ("cpu_cores", 0, False, "below minimum"),
            ("cpu_cores", 64, False, "exceeds maximum"),
            ("memory_mb", 2000, True, ""),
            ("memory_mb", 512, True, ""),
            ("memory_mb", 32768, True, ""),
            ("memory_mb", 511, False, "below minimum"),
            ("memory_mb", 32769, False, "exceeds maximum"),
            ("charges", -10, True, ""),
            ("charges", 0, True, ""),
            ("charges", 10, True, ""),
            ("charges", -11, False, "below minimum"),
            ("charges", 11, False, "exceeds maximum"),
            ("spin", 0, True, ""),
            ("spin", 10, True, ""),
            ("spin", -1, False, "below minimum"),
            ("spin", 11, False, "exceeds maximum"),
            ("unknown_param", 999, True, ""),
        ],
    )
    def test_validate_parameter_value(
        self,
        param,
        value,
        expected_valid,
        expected_error,
    ):
        """Test parameter value bounds and unknown-parameter behavior."""
        is_valid, error = validate_parameter_value(param, value)
        assert is_valid is expected_valid
        if expected_valid:
            assert error == ""
        else:
            assert expected_error in error


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
            assert defaults['memory_mb'] >= 512
            assert defaults['memory_mb'] <= 32768

    def test_constraint_min_less_than_max(self):
        """Test that min is always less than max in constraints."""
        for param_name, constraint in PARAMETER_CONSTRAINTS.items():
            if 'min' in constraint and 'max' in constraint:
                assert constraint['min'] < constraint['max'], \
                    f"Parameter {param_name}: min must be less than max"


class TestValidateParametersForMethod:
    """Tests for validate_parameters_for_method function."""

    def test_valid_dft_parameters(self):
        """Test that valid DFT parameters pass validation."""
        params = {
            'xyz': 'H 0 0 0\nH 0 0 0.74',
            'calculation_method': 'DFT',
            'basis_function': '6-31G(d)',
            'exchange_correlation': 'B3LYP',
            'charges': 0,
            'spin': 0,
            'optimize_geometry': True
        }

        is_valid, error = validate_parameters_for_method('DFT', params)
        assert is_valid is True
        assert error == ''

    def test_dft_rejects_casci_parameters(self):
        """Test that DFT rejects CASCI-specific parameters (ncas, nelecas)."""
        params = {
            'xyz': 'H 0 0 0\nH 0 0 0.74',
            'calculation_method': 'DFT',
            'basis_function': '6-31G(d)',
            'ncas': 4,  # Not applicable to DFT
            'nelecas': 4  # Not applicable to DFT
        }

        is_valid, error = validate_parameters_for_method('DFT', params)
        assert is_valid is False
        assert 'ncas' in error
        assert 'not applicable' in error.lower()
        assert 'CASCI' in error or 'CASSCF' in error

    def test_dft_rejects_tddft_parameters(self):
        """Test that DFT rejects TDDFT-specific parameters."""
        params = {
            'xyz': 'H 0 0 0\nH 0 0 0.74',
            'calculation_method': 'DFT',
            'tddft_nstates': 10,  # Not applicable to DFT
            'tddft_analyze_nto': True  # Not applicable to DFT
        }

        is_valid, error = validate_parameters_for_method('DFT', params)
        assert is_valid is False
        assert 'tddft_nstates' in error or 'tddft_analyze_nto' in error
        assert 'not applicable' in error.lower()

    def test_valid_casci_parameters(self):
        """Test that valid CASCI parameters pass validation."""
        params = {
            'xyz': 'H 0 0 0\nH 0 0 0.74',
            'calculation_method': 'CASCI',
            'basis_function': '6-31G(d)',
            'ncas': 4,
            'nelecas': 4,
            'natorb': True,
            'max_cycle_micro': 3
        }

        is_valid, error = validate_parameters_for_method('CASCI', params)
        assert is_valid is True
        assert error == ''

    def test_casci_rejects_out_of_range_active_space(self):
        """Test that CASCI validation enforces active space bounds."""
        params = {
            'xyz': 'H 0 0 0\nH 0 0 0.74',
            'calculation_method': 'CASCI',
            'ncas': 21,
            'nelecas': 41
        }

        is_valid, error = validate_parameters_for_method('CASCI', params)
        assert is_valid is False
        assert 'ncas' in error
        assert 'nelecas' in error
        assert 'exceeds maximum' in error

    def test_universal_parameters_still_validate_bounds(self):
        """Test that universal parameters are accepted only within their value bounds."""
        params = {
            'xyz': 'H 0 0 0\nH 0 0 0.74',
            'calculation_method': 'DFT',
            'basis_function': '6-31G(d)',
            'cpu_cores': 64,
            'memory_mb': 64,
            'charges': 11,
            'spin': 11
        }

        is_valid, error = validate_parameters_for_method('DFT', params)
        assert is_valid is False
        assert 'cpu_cores' in error
        assert 'memory_mb' in error
        assert 'charges' in error
        assert 'spin' in error

    def test_casci_rejects_tddft_parameters(self):
        """Test that CASCI rejects TDDFT-specific parameters."""
        params = {
            'xyz': 'H 0 0 0\nH 0 0 0.74',
            'calculation_method': 'CASCI',
            'ncas': 4,
            'nelecas': 4,
            'tddft_nstates': 10  # Not applicable to CASCI
        }

        is_valid, error = validate_parameters_for_method('CASCI', params)
        assert is_valid is False
        assert 'tddft_nstates' in error
        assert 'not applicable' in error.lower()

    def test_valid_tddft_parameters(self):
        """Test that valid TDDFT parameters pass validation."""
        params = {
            'xyz': 'H 0 0 0\nH 0 0 0.74',
            'calculation_method': 'TDDFT',
            'basis_function': '6-31G(d)',
            'exchange_correlation': 'B3LYP',
            'tddft_nstates': 10,
            'tddft_method': 'TDDFT',
            'tddft_analyze_nto': False
        }

        is_valid, error = validate_parameters_for_method('TDDFT', params)
        assert is_valid is True
        assert error == ''

    def test_tddft_rejects_optimize_geometry(self):
        """Test that TDDFT rejects optimize_geometry parameter."""
        params = {
            'xyz': 'H 0 0 0\\nH 0 0 0.74',
            'calculation_method': 'TDDFT',
            'exchange_correlation': 'B3LYP',
            'tddft_nstates': 10,
            'optimize_geometry': True  # Not applicable for TDDFT
        }

        is_valid, error = validate_parameters_for_method('TDDFT', params)
        assert is_valid is False
        assert 'optimize_geometry' in error
        assert 'not applicable' in error.lower()



    def test_ccsd_accepts_frozen_core(self):
        """Test that CCSD accepts frozen_core parameter."""
        params = {
            'xyz': 'H 0 0 0\nH 0 0 0.74',
            'calculation_method': 'CCSD',
            'basis_function': 'cc-pVDZ',
            'frozen_core': True
        }

        is_valid, error = validate_parameters_for_method('CCSD', params)
        assert is_valid is True
        assert error == ''

    def test_dft_rejects_frozen_core(self):
        """Test that DFT rejects frozen_core parameter."""
        params = {
            'xyz': 'H 0 0 0\nH 0 0 0.74',
            'calculation_method': 'DFT',
            'exchange_correlation': 'B3LYP',
            'frozen_core': True  # Not applicable to DFT
        }

        is_valid, error = validate_parameters_for_method('DFT', params)
        assert is_valid is False
        assert 'frozen_core' in error
        assert 'not applicable' in error.lower()

    def test_universal_parameters_accepted_by_all_methods(self):
        """Test that universal parameters are accepted by all methods."""
        universal_params = {
            'xyz': 'H 0 0 0\nH 0 0 0.74',
            'calculation_method': 'DFT',
            'basis_function': '6-31G(d)',
            'charges': -1,
            'spin': 1,
            'solvent_method': 'ief-pcm',
            'solvent': 'water',
            'name': 'test molecule',
            'cpu_cores': 4,
            'memory_mb': 2000
        }

        # Test each method accepts universal parameters
        for method in ['DFT', 'HF', 'MP2', 'CCSD', 'TDDFT', 'CASCI', 'CASSCF']:
            test_params = {**universal_params, 'calculation_method': method}
            is_valid, error = validate_parameters_for_method(method, test_params)
            assert is_valid is True, f"Method {method} should accept universal parameters. Error: {error}"

    def test_none_values_are_ignored(self):
        """Test that None values (unprovided parameters) are ignored."""
        params = {
            'xyz': 'H 0 0 0\nH 0 0 0.74',
            'calculation_method': 'DFT',
            'basis_function': '6-31G(d)',
            'ncas': None,  # None should be ignored
            'nelecas': None,  # None should be ignored
            'tddft_nstates': None  # None should be ignored
        }

        is_valid, error = validate_parameters_for_method('DFT', params)
        assert is_valid is True
        assert error == ''

    def test_multiple_invalid_parameters_in_error_message(self):
        """Test that error message includes all invalid parameters."""
        params = {
            'xyz': 'H 0 0 0\nH 0 0 0.74',
            'calculation_method': 'DFT',
            'ncas': 4,  # Invalid
            'nelecas': 4,  # Invalid
            'tddft_nstates': 10  # Invalid
        }

        is_valid, error = validate_parameters_for_method('DFT', params)
        assert is_valid is False
        # All three invalid parameters should be mentioned
        assert 'ncas' in error
        assert 'nelecas' in error or 'tddft_nstates' in error

    def test_hf_rejects_exchange_correlation(self):
        """Test that HF rejects exchange_correlation parameter."""
        params = {
            'xyz': 'H 0 0 0\nH 0 0 0.74',
            'calculation_method': 'HF',
            'basis_function': '6-31G(d)',
            'exchange_correlation': 'B3LYP'  # Not applicable to HF
        }

        is_valid, error = validate_parameters_for_method('HF', params)
        assert is_valid is False
        assert 'exchange_correlation' in error
        assert 'not applicable' in error.lower()
        assert 'DFT' in error or 'TDDFT' in error

    def test_valid_hf_parameters(self):
        """Test that valid HF parameters pass validation."""
        params = {
            'xyz': 'H 0 0 0\nH 0 0 0.74',
            'calculation_method': 'HF',
            'basis_function': '6-31G(d)',
            'charges': 0,
            'spin': 0,
            'optimize_geometry': True
        }

        is_valid, error = validate_parameters_for_method('HF', params)
        assert is_valid is True
        assert error == ''

    def test_mp2_rejects_exchange_correlation(self):
        """Test that MP2 rejects exchange_correlation parameter."""
        params = {
            'xyz': 'H 0 0 0\nH 0 0 0.74',
            'calculation_method': 'MP2',
            'basis_function': '6-31G(d)',
            'exchange_correlation': 'PBE'  # Not applicable to MP2
        }

        is_valid, error = validate_parameters_for_method('MP2', params)
        assert is_valid is False
        assert 'exchange_correlation' in error
        assert 'not applicable' in error.lower()
