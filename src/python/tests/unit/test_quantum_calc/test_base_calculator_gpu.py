"""
Unit tests for BaseCalculator GPU detection helpers.
"""

from types import SimpleNamespace
from unittest.mock import MagicMock

import pytest
import quantum_calc.base_calculator as base_calculator
from quantum_calc.base_calculator import BaseCalculator
from quantum_calc.exceptions import CalculationError


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


class FailingKernel:
    """Mean-field stand-in that fails during kernel execution."""

    def kernel(self):
        raise RuntimeError("CUDA execution failed")


def test_gpu_disabled_returns_false(monkeypatch):
    """
    GIVEN GPU acceleration disabled in settings
    WHEN _is_gpu4pyscf_available is called
    THEN it returns False without probing modules
    """
    calculator = DummyCalculator(optimize_geometry=False)
    monkeypatch.setattr(calculator, "_is_gpu_acceleration_enabled", lambda: False)

    find_spec = MagicMock()
    monkeypatch.setattr(base_calculator.util, "find_spec", find_spec)

    assert calculator._is_gpu4pyscf_available() is False
    find_spec.assert_not_called()


def test_gpu_unavailable_on_non_linux(monkeypatch):
    """
    GIVEN GPU acceleration enabled on non-Linux
    WHEN _is_gpu4pyscf_available is called
    THEN it returns False
    """
    calculator = DummyCalculator(optimize_geometry=False)
    monkeypatch.setattr(calculator, "_is_gpu_acceleration_enabled", lambda: True)
    monkeypatch.setattr(base_calculator.sys, "platform", "darwin")

    find_spec = MagicMock()
    monkeypatch.setattr(base_calculator.util, "find_spec", find_spec)

    assert calculator._is_gpu4pyscf_available() is False
    find_spec.assert_not_called()


def test_gpu_requires_supported_cuda(monkeypatch):
    """
    GIVEN Linux with GPU acceleration enabled but CUDA not detected
    WHEN _is_gpu4pyscf_available is called
    THEN it returns False without module probing
    """
    calculator = DummyCalculator(optimize_geometry=False)
    monkeypatch.setattr(calculator, "_is_gpu_acceleration_enabled", lambda: True)
    monkeypatch.setattr(base_calculator.sys, "platform", "linux")
    monkeypatch.setattr(base_calculator.shutil, "which", lambda _: None)

    find_spec = MagicMock()
    monkeypatch.setattr(base_calculator.util, "find_spec", find_spec)

    assert calculator._is_gpu4pyscf_available() is False
    find_spec.assert_not_called()


def test_gpu_available_when_cuda_supported_and_module_present(monkeypatch):
    """
    GIVEN Linux with supported CUDA and gpu4pyscf installed
    WHEN _is_gpu4pyscf_available is called
    THEN it returns True
    """
    calculator = DummyCalculator(optimize_geometry=False)
    monkeypatch.setattr(calculator, "_is_gpu_acceleration_enabled", lambda: True)
    monkeypatch.setattr(base_calculator.sys, "platform", "linux")
    monkeypatch.setattr(base_calculator.shutil, "which", lambda _: "/usr/local/cuda/bin/nvcc")

    nvcc_output = SimpleNamespace(
        returncode=0,
        stdout="Cuda compilation tools, release 12.4, V12.4.0",
        stderr="",
    )
    monkeypatch.setattr(base_calculator.subprocess, "run", lambda *args, **kwargs: nvcc_output)
    monkeypatch.setattr(base_calculator.util, "find_spec", lambda _: object())

    assert calculator._is_gpu4pyscf_available() is True


def test_gpu_required_raises_when_enabled_but_unavailable(monkeypatch):
    """
    GIVEN GPU acceleration enabled on an unsupported platform
    WHEN GPU execution is required
    THEN a calculation error is raised instead of falling back to CPU
    """
    calculator = DummyCalculator(optimize_geometry=False)
    monkeypatch.setattr(calculator, "_is_gpu_acceleration_enabled", lambda: True)
    monkeypatch.setattr(base_calculator.sys, "platform", "darwin")

    with pytest.raises(CalculationError, match="GPU acceleration is enabled"):
        calculator._require_gpu4pyscf_available()


def test_gpu_runtime_failure_is_reported_as_calculation_error():
    """
    GIVEN a GPU-backed mean-field object
    WHEN the PySCF kernel fails during execution
    THEN the error is surfaced as a GPU calculation error
    """
    calculator = DummyCalculator(optimize_geometry=False)
    calculator.mf = FailingKernel()
    calculator.gpu_enabled = True

    with pytest.raises(CalculationError, match="GPU4PySCF calculation failed"):
        calculator._run_base_scf_calculation()
