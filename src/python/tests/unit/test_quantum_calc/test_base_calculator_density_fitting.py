"""
Unit tests for BaseCalculator density fitting helpers.
"""

from types import SimpleNamespace

from quantum_calc.base_calculator import BaseCalculator


class DummyCalculator(BaseCalculator):
    """Minimal concrete calculator for testing BaseCalculator helpers."""

    def _perform_specific_calculation(self, base_energy: float):
        return {}

    def _create_scf_method(self, mol):
        return None

    def _apply_solvent_effects(self, mf):
        return mf

    def _get_base_method_description(self) -> str:
        return "Dummy"


def test_resolve_actual_auxiliary_basis_returns_none_when_density_fitting_disabled():
    calculator = DummyCalculator(optimize_geometry=False)
    calculator.density_fitting = False
    calculator.mf = SimpleNamespace(with_df=SimpleNamespace(auxbasis='weigend'))

    assert calculator._resolve_actual_auxiliary_basis() is None


def test_resolve_actual_auxiliary_basis_prefers_explicit_auxbasis_name():
    calculator = DummyCalculator(optimize_geometry=False)
    calculator.density_fitting = True
    calculator.mf = SimpleNamespace(with_df=SimpleNamespace(auxbasis='def2-tzvp-jkfit'))

    assert calculator._resolve_actual_auxiliary_basis() == 'def2-tzvp-jkfit'


def test_resolve_actual_auxiliary_basis_collapses_uniform_element_mapping():
    calculator = DummyCalculator(optimize_geometry=False)
    calculator.density_fitting = True
    calculator.mf = SimpleNamespace(
        with_df=SimpleNamespace(
            auxbasis=None,
            auxmol=SimpleNamespace(
                basis={'O': 'cc-pvdz-jkfit', 'H': 'cc-pvdz-jkfit'}
            ),
        )
    )

    assert calculator._resolve_actual_auxiliary_basis() == 'cc-pvdz-jkfit'


def test_resolve_actual_auxiliary_basis_formats_mixed_element_mapping():
    calculator = DummyCalculator(optimize_geometry=False)
    calculator.density_fitting = True
    calculator.mf = SimpleNamespace(
        with_df=SimpleNamespace(
            auxbasis=None,
            auxmol=SimpleNamespace(
                basis={'C': 'def2-svp-jkfit', 'H': 'cc-pvdz-jkfit'}
            ),
        )
    )

    assert (
        calculator._resolve_actual_auxiliary_basis()
        == 'C: def2-svp-jkfit, H: cc-pvdz-jkfit'
    )
