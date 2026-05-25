"""Unit tests for checkpoint resume callbacks."""

from __future__ import annotations

import sys
from types import ModuleType, SimpleNamespace
from unittest.mock import Mock

from quantum_calc._checkpoint_resume import CheckpointResumeMixin
from quantum_calc import _worker_runtime


class DummyCheckpointResumeCalculator(CheckpointResumeMixin):
    """Minimal calculator exposing the checkpoint resume mixin callback."""

    def __init__(self, working_dir: str) -> None:
        self.working_dir = working_dir
        self.file_manager = Mock()


class DummyMol:
    """Small molecule stub with the PySCF geometry API used by the callback."""

    natm = 2

    def atom_symbol(self, index: int) -> str:
        return ["H", "O"][index]

    def atom_coords(self, unit: str = "ANG") -> list[list[float]]:
        assert unit == "ANG"
        return [
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.9572],
        ]


def test_geometry_optimization_callback_with_dict_envs_saves_trajectory_step(
    tmp_path,
):
    """
    GIVEN PySCF provides geometry optimizer envs as a dict
    WHEN the geometry optimization callback runs
    THEN the current geometry is saved as an XYZ trajectory step
    """
    calculator = DummyCheckpointResumeCalculator(str(tmp_path))
    mock_mol = DummyMol()
    envs = {"mol": mock_mol, "cycle": 3}

    result = calculator._geometry_optimization_callback(envs)

    assert result is False
    save_step = calculator.file_manager.save_geometry_trajectory_step
    save_step.assert_called_once()
    working_dir, step_num, geometry_xyz = save_step.call_args.args
    assert working_dir == str(tmp_path)
    assert step_num == 3
    assert geometry_xyz.splitlines()[0] == "2"
    assert geometry_xyz.count("Optimization step") == 1
    assert "Optimization step 3" in geometry_xyz


class FakeThreadPoolLimits:
    """No-op context manager for worker runtime unit tests."""

    def __init__(self, limits: int) -> None:
        self.limits = limits

    def __enter__(self) -> None:
        return None

    def __exit__(
        self,
        exc_type: object,
        exc_value: object,
        traceback: object,
    ) -> bool:
        return False


class FakePySCFLib:
    """Small PySCF lib.num_threads stub."""

    def __init__(self) -> None:
        self.thread_count = 1

    def num_threads(self, value: int | None = None) -> int:
        if value is not None:
            self.thread_count = value
        return self.thread_count


class FakeGeometryResumeCalculator:
    """Calculator stub that records setup inputs."""

    def __init__(self, last_geometry: str) -> None:
        self.file_manager = Mock()
        self.file_manager.load_last_geometry.return_value = last_geometry
        self.parsed_xyz_values: list[str] = []
        self.setup_atoms = None
        self.resume_from_checkpoint_called = False

    def parse_xyz(self, xyz: str) -> list[tuple[str, str]]:
        self.parsed_xyz_values.append(xyz)
        return [("parsed_from", xyz)]

    def setup_calculation(self, atoms: list[tuple[str, str]], **kwargs: object) -> None:
        self.setup_atoms = atoms

    def resume_from_checkpoint(self) -> None:
        self.resume_from_checkpoint_called = True

    def run_calculation(self) -> dict:
        return {"energy": -1.0}


def test_calculation_worker_geometry_resume_uses_last_geometry_for_setup(
    tmp_path,
    monkeypatch,
):
    """
    GIVEN a paused geometry optimization has a saved trajectory geometry
    WHEN the worker resumes the calculation
    THEN setup_calculation receives atoms parsed from the last trajectory geometry
    """
    original_xyz = "2\noriginal\nH 0 0 0\nH 0 0 1"
    last_geometry = "2\nlast\nH 0 0 0\nH 0 0 2"
    pause_state = {"calculation_phase": "geometry_optimization"}
    fake_calculator = FakeGeometryResumeCalculator(last_geometry)
    fake_pyscf = ModuleType("pyscf")
    fake_pyscf.lib = FakePySCFLib()
    fake_threadpoolctl = ModuleType("threadpoolctl")
    fake_threadpoolctl.threadpool_info = lambda: []
    fake_threadpoolctl.threadpool_limits = FakeThreadPoolLimits

    calc_dir = tmp_path / "calc-geometry-resume"
    calc_dir.mkdir()

    monkeypatch.setitem(sys.modules, "pyscf", fake_pyscf)
    monkeypatch.setitem(sys.modules, "threadpoolctl", fake_threadpoolctl)
    monkeypatch.setattr(
        "quantum_calc.get_current_settings",
        lambda: SimpleNamespace(calculations_directory=str(tmp_path)),
    )
    monkeypatch.setattr(_worker_runtime, "_setup_worker_environment", lambda *_: (1, 512))
    monkeypatch.setattr(_worker_runtime, "_check_casci_dependencies", lambda *_: None)
    monkeypatch.setattr(
        _worker_runtime,
        "_import_calculator_classes",
        lambda *_: {"DFT": object},
    )
    monkeypatch.setattr(
        _worker_runtime,
        "_create_calculator_instance",
        lambda *args: fake_calculator,
    )

    success, error = _worker_runtime.calculation_worker(
        "calc-geometry-resume",
        {
            "name": "Water",
            "xyz": original_xyz,
            "calculation_method": "DFT",
            "basis_function": "sto-3g",
            "charges": 0,
            "spin": 0,
            "solvent_method": None,
            "solvent": None,
            "exchange_correlation": "b3lyp",
            "resume_from_pause": True,
            "pause_state": pause_state,
        },
    )

    assert success is True
    assert error is None
    assert fake_calculator.parsed_xyz_values == [last_geometry]
    assert fake_calculator.setup_atoms == [("parsed_from", last_geometry)]
    assert fake_calculator.resume_from_checkpoint_called is True
