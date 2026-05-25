"""
Unit tests for BaseCalculator GPU detection helpers.
"""

import sys
import types
from types import SimpleNamespace
from unittest.mock import MagicMock

import pytest
from pyscf import gto
import quantum_calc.base_calculator as base_calculator
from quantum_calc.base_calculator import BaseCalculator
from quantum_calc.dft_calculator import DFTCalculator
from quantum_calc.exceptions import CalculationError, PauseRequestedException
from quantum_calc.hf_calculator import HFCalculator
from quantum_calc.tddft_calculator import TDDFTCalculator


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
    THEN it raises instead of falling back to CPU
    """
    calculator = DummyCalculator(optimize_geometry=False)
    monkeypatch.setattr(calculator, "_is_gpu_acceleration_enabled", lambda: True)
    monkeypatch.setattr(base_calculator.sys, "platform", "darwin")

    with pytest.raises(CalculationError, match="GPU acceleration is enabled"):
        calculator._require_gpu4pyscf_available()


def test_hf_gpu_setup_failure_raises(tmp_path, monkeypatch):
    """
    GIVEN GPU prerequisites are available but GPU4PySCF HF setup fails
    WHEN the SCF method is created
    THEN the setup error is reported instead of falling back to CPU
    """
    gpu_scf = types.ModuleType("gpu4pyscf.scf")
    gpu_scf.RHF = MagicMock(side_effect=RuntimeError("GPU setup failed"))
    gpu_scf.UHF = MagicMock(side_effect=RuntimeError("GPU setup failed"))
    gpu4pyscf = types.ModuleType("gpu4pyscf")
    gpu4pyscf.scf = gpu_scf
    monkeypatch.setitem(sys.modules, "gpu4pyscf", gpu4pyscf)
    monkeypatch.setitem(sys.modules, "gpu4pyscf.scf", gpu_scf)

    calculator = HFCalculator(working_dir=str(tmp_path), optimize_geometry=False)
    calculator.results["spin"] = 0
    monkeypatch.setattr(calculator, "_require_gpu4pyscf_available", lambda: True)
    mol = gto.M(atom="H 0 0 0; H 0 0 0.74", basis="sto-3g", verbose=0)

    with pytest.raises(CalculationError, match="GPU4PySCF HF setup failed"):
        calculator._create_scf_method(mol)


def test_dft_gpu_setup_failure_raises(tmp_path, monkeypatch):
    """
    GIVEN GPU prerequisites are available but GPU4PySCF DFT setup fails
    WHEN the SCF method is created
    THEN the setup error is reported instead of falling back to CPU
    """
    gpu_dft = types.ModuleType("gpu4pyscf.dft")
    gpu_dft.RKS = MagicMock(side_effect=RuntimeError("GPU setup failed"))
    gpu_dft.UKS = MagicMock(side_effect=RuntimeError("GPU setup failed"))
    gpu4pyscf = types.ModuleType("gpu4pyscf")
    gpu4pyscf.dft = gpu_dft
    monkeypatch.setitem(sys.modules, "gpu4pyscf", gpu4pyscf)
    monkeypatch.setitem(sys.modules, "gpu4pyscf.dft", gpu_dft)

    calculator = DFTCalculator(working_dir=str(tmp_path), optimize_geometry=False)
    calculator.results["spin"] = 0
    calculator.xc_functional = "B3LYP"
    monkeypatch.setattr(calculator, "_require_gpu4pyscf_available", lambda: True)
    mol = gto.M(atom="H 0 0 0; H 0 0 0.74", basis="sto-3g", verbose=0)

    with pytest.raises(CalculationError, match="GPU4PySCF DFT setup failed"):
        calculator._create_scf_method(mol)


def test_gpu_base_scf_failure_raises_without_cpu_fallback(monkeypatch):
    """
    GIVEN GPU execution is active and base SCF kernel fails
    WHEN _run_base_scf_calculation is called
    THEN the GPU failure is reported without CPU fallback
    """
    calculator = DummyCalculator(optimize_geometry=False)
    calculator.gpu_enabled = True
    calculator.mf = SimpleNamespace(kernel=MagicMock(side_effect=RuntimeError("gpu failed")))

    with pytest.raises(CalculationError, match="GPU4PySCF calculation failed"):
        calculator._run_base_scf_calculation()


def test_gpu_base_scf_pause_propagates_without_cpu_fallback(monkeypatch):
    """
    GIVEN GPU execution is active and base SCF kernel requests pause
    WHEN _run_base_scf_calculation is called
    THEN the pause exception propagates without CPU fallback
    """
    calculator = DummyCalculator(optimize_geometry=False)
    pause_error = PauseRequestedException("pause requested")
    calculator.gpu_enabled = True
    calculator.mf = SimpleNamespace(kernel=MagicMock(side_effect=pause_error))

    with pytest.raises(PauseRequestedException) as exc_info:
        calculator._run_base_scf_calculation()

    assert exc_info.value is pause_error


def test_tddft_gpu_kernel_pause_propagates_without_cpu_fallback(tmp_path, monkeypatch):
    """
    GIVEN GPU TDDFT execution is active and TDDFT kernel requests pause
    WHEN _perform_specific_calculation is called
    THEN the pause exception propagates without CPU fallback
    """
    calculator = TDDFTCalculator(working_dir=str(tmp_path), optimize_geometry=False)
    pause_error = PauseRequestedException("pause requested")
    mytd = SimpleNamespace(
        nstates=1,
        kernel=MagicMock(side_effect=pause_error),
    )
    calculator.gpu_enabled = True
    calculator.tddft_nstates = 1
    calculator.tddft_method = "TDDFT"
    calculator.mf = SimpleNamespace(
        mo_energy=[-0.5, -0.1, 0.2, 0.4],
        mo_occ=[2, 2, 0, 0],
        TDDFT=MagicMock(return_value=mytd),
    )
    with pytest.raises(PauseRequestedException) as exc_info:
        calculator._perform_specific_calculation(-1.0)

    assert exc_info.value is pause_error


def test_tddft_gpu_kernel_failure_raises_without_cpu_fallback(tmp_path):
    """
    GIVEN GPU TDDFT execution is active and TDDFT kernel fails
    WHEN _perform_specific_calculation is called
    THEN the GPU failure is reported without CPU fallback
    """
    calculator = TDDFTCalculator(working_dir=str(tmp_path), optimize_geometry=False)
    mytd = SimpleNamespace(
        nstates=1,
        kernel=MagicMock(side_effect=RuntimeError("td gpu failed")),
    )
    calculator.gpu_enabled = True
    calculator.tddft_nstates = 1
    calculator.tddft_method = "TDDFT"
    calculator.mf = SimpleNamespace(
        mo_energy=[-0.5, -0.1, 0.2, 0.4],
        mo_occ=[2, 2, 0, 0],
        TDDFT=MagicMock(return_value=mytd),
    )

    with pytest.raises(CalculationError, match="GPU4PySCF TDDFT calculation failed"):
        calculator._perform_specific_calculation(-1.0)
