"""Unit tests for CalculationRepository status validation."""

import json
from datetime import datetime
from pathlib import Path

import quantum_calc._calculation_repository as calculation_repository
from quantum_calc._calculation_repository import CalculationRepository


def test_read_calculation_status_invalid_value_returns_error(tmp_path):
    """
    GIVEN status.json contains an invalid status value
    WHEN read_calculation_status is called
    THEN it should return "error"
    """
    manager = CalculationRepository(base_dir=str(tmp_path))
    calc_dir = tmp_path / "calc_invalid_status"
    calc_dir.mkdir()

    with open(calc_dir / "status.json", "w") as f:
        json.dump({"status": "invalid_status"}, f)

    assert manager.read_calculation_status(str(calc_dir)) == "error"


def test_read_calculation_status_details_invalid_value_returns_error(tmp_path):
    """
    GIVEN status.json contains an invalid status value
    WHEN read_calculation_status_details is called
    THEN it should return ("error", None)
    """
    manager = CalculationRepository(base_dir=str(tmp_path))
    calc_dir = tmp_path / "calc_invalid_status_details"
    calc_dir.mkdir()

    with open(calc_dir / "status.json", "w") as f:
        json.dump({"status": "invalid_status"}, f)

    assert manager.read_calculation_status_details(str(calc_dir)) == ("error", None)


def test_read_calculation_status_valid_value_returns_as_is(tmp_path):
    """
    GIVEN status.json contains a valid status value
    WHEN status readers are called
    THEN the status should be returned as-is
    """
    manager = CalculationRepository(base_dir=str(tmp_path))
    calc_dir = tmp_path / "calc_valid_status"
    calc_dir.mkdir()

    with open(calc_dir / "status.json", "w") as f:
        json.dump({"status": "running"}, f)

    assert manager.read_calculation_status(str(calc_dir)) == "running"
    assert manager.read_calculation_status_details(str(calc_dir)) == ("running", None)


def test_create_calculation_dir_same_name_same_second_creates_unique_directories(
    tmp_path, monkeypatch
):
    """
    GIVEN two calculations with the same molecule name are created in the same second
    WHEN create_calculation_dir is called repeatedly
    THEN each call should create a distinct directory
    """
    manager = CalculationRepository(base_dir=str(tmp_path))

    class FixedDatetime:
        @classmethod
        def now(cls):
            return datetime(2026, 5, 20, 12, 34, 56, 123456)

    monkeypatch.setattr(calculation_repository, "datetime", FixedDatetime)

    first_dir = manager.create_calculation_dir("Water")
    second_dir = manager.create_calculation_dir("Water")

    assert first_dir != second_dir
    assert (tmp_path / Path(first_dir).name).is_dir()
    assert (tmp_path / Path(second_dir).name).is_dir()
