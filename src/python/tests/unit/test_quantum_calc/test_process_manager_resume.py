"""
Unit tests for CalculationProcessManager pause/resume persistence.
"""

from types import SimpleNamespace

from quantum_calc._calculation_repository import CalculationRepository
from quantum_calc.process_manager import CalculationProcessManager


def test_resume_calculation_persists_waiting_status(tmp_path, monkeypatch):
    """
    GIVEN a paused calculation is resumed but queued
    WHEN resume_calculation receives a waiting submit result
    THEN status.json is updated to waiting with the queue reason
    """
    repository = CalculationRepository(base_dir=str(tmp_path))
    calc_id = "paused-calc"
    calc_dir = tmp_path / calc_id
    calc_dir.mkdir()
    repository.save_calculation_parameters(
        str(calc_dir),
        {
            "name": "Paused Calc",
            "created_at": "2026-05-20T00:00:00",
            "calculation_method": "HF",
        },
    )
    repository.save_calculation_status(str(calc_dir), "paused")

    manager = object.__new__(CalculationProcessManager)
    manager._shutdown = True
    monkeypatch.setattr(
        "quantum_calc.get_current_settings",
        lambda: SimpleNamespace(calculations_directory=str(tmp_path)),
    )
    monkeypatch.setattr(
        manager,
        "submit_calculation",
        lambda calculation_id, params: (True, "waiting", "All slots are busy"),
    )

    result = manager.resume_calculation(calc_id)

    assert result == {
        "calculation_id": calc_id,
        "status": "waiting",
        "waiting_reason": "All slots are busy",
    }
    assert repository.read_calculation_status_details(str(calc_dir)) == (
        "waiting",
        "All slots are busy",
    )
