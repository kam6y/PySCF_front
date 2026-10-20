"""Unit tests for CalculationRepository status validation."""

import json
import os
from datetime import datetime
from pathlib import Path

import pytest

import quantum_calc._calculation_repository as calculation_repository
from quantum_calc._calculation_repository import CalculationRepository


@pytest.mark.parametrize(
    (
        "method_name",
        "filename",
        "existing_data",
        "new_args",
        "expected_existing_data",
    ),
    [
        (
            "save_calculation_parameters",
            "parameters.json",
            {"name": "old calculation", "basis_function": "sto-3g"},
            ({"name": "new calculation", "basis_function": "6-31g"},),
            {"name": "old calculation", "basis_function": "sto-3g"},
        ),
        (
            "save_calculation_results",
            "results.json",
            {"energy": -75.0, "success": True},
            ({"energy": -74.0, "success": False},),
            {"energy": -75.0, "success": True},
        ),
        (
            "save_calculation_status",
            "status.json",
            {"status": "completed", "updated_at": "2026-05-20T00:00:00"},
            ("running",),
            {"status": "completed", "updated_at": "2026-05-20T00:00:00"},
        ),
    ],
)
def test_save_major_json_preserves_existing_file_when_replace_fails(
    tmp_path,
    monkeypatch,
    method_name,
    filename,
    existing_data,
    new_args,
    expected_existing_data,
):
    """
    GIVEN a major JSON file already exists
    WHEN atomic replacement fails while saving new content
    THEN the existing JSON file should remain unchanged
    """
    manager = CalculationRepository(base_dir=str(tmp_path))
    calc_dir = tmp_path / "calc_atomic_failure"
    calc_dir.mkdir()
    target_file = calc_dir / filename
    target_file.write_text(json.dumps(existing_data), encoding="utf-8")

    def fail_replace(src, dst):
        raise OSError("simulated replace failure")

    monkeypatch.setattr(os, "replace", fail_replace)

    with pytest.raises(OSError, match="simulated replace failure"):
        getattr(manager, method_name)(str(calc_dir), *new_args)

    assert json.loads(target_file.read_text(encoding="utf-8")) == expected_existing_data


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


def test_read_calculation_status_corrupt_json_returns_error(tmp_path):
    """
    GIVEN status.json exists but is not valid JSON
    WHEN status readers are called
    THEN they should treat the calculation as error instead of pending
    """
    manager = CalculationRepository(base_dir=str(tmp_path))
    calc_dir = tmp_path / "calc_corrupt_status"
    calc_dir.mkdir()

    (calc_dir / "status.json").write_text("{not valid json", encoding="utf-8")

    assert manager.read_calculation_status(str(calc_dir)) == "error"
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


def test_resolve_calculation_path_valid_id_returns_path_under_base_dir(tmp_path):
    """
    GIVEN a valid calculation ID
    WHEN resolve_calculation_path is called
    THEN it should return the path under the repository base directory
    """
    manager = CalculationRepository(base_dir=str(tmp_path))

    resolved_path = manager.resolve_calculation_path("calc-123")

    assert resolved_path == tmp_path.resolve() / "calc-123"


@pytest.mark.parametrize("calculation_id", ["", ".", "..", "a/b", r"a\b"])
def test_resolve_calculation_path_rejects_unsafe_ids(tmp_path, calculation_id):
    """
    GIVEN an unsafe calculation ID
    WHEN resolve_calculation_path is called
    THEN it should reject the ID before filesystem access
    """
    manager = CalculationRepository(base_dir=str(tmp_path))

    with pytest.raises(ValueError, match="Invalid calculation ID"):
        manager.resolve_calculation_path(calculation_id)
