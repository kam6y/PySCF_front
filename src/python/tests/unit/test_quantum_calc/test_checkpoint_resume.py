"""Unit tests for checkpoint resume callbacks."""

from __future__ import annotations

from unittest.mock import Mock

from quantum_calc._checkpoint_resume import CheckpointResumeMixin


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
