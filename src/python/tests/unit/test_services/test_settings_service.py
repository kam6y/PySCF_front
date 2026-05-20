from pathlib import Path
from types import SimpleNamespace
from typing import Any
from unittest.mock import Mock

import pytest

import services.settings_service as settings_service_module
from services.exceptions import ValidationError
from services.settings_service import SettingsService


class _UpdatedSettings(SimpleNamespace):
    def model_dump(self, mode: str | None = None) -> dict[str, Any]:
        return self.__dict__


def test_update_settings_with_directory_migration_conflict_keeps_current_settings(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """
    GIVEN the new calculations directory already contains a conflicting calculation
    WHEN settings update requests a calculations_directory change
    THEN settings and QuantumService are not switched to the new directory
    """
    old_dir = tmp_path / "old_calculations"
    new_dir = tmp_path / "new_calculations"
    old_calc_a = old_dir / "calc_a"
    old_calc_b = old_dir / "calc_b"
    existing_new_calc_b = new_dir / "calc_b"

    old_calc_a.mkdir(parents=True)
    old_calc_b.mkdir()
    existing_new_calc_b.mkdir(parents=True)
    (old_calc_a / "metadata.json").write_text("old a", encoding="utf-8")
    (old_calc_b / "metadata.json").write_text("old b", encoding="utf-8")
    (existing_new_calc_b / "metadata.json").write_text("new b", encoding="utf-8")

    current_settings = SimpleNamespace(calculations_directory=str(old_dir))
    updated_settings = _UpdatedSettings(
        calculations_directory=str(new_dir),
        max_parallel_instances=2,
        max_cpu_utilization_percent=90.0,
        max_memory_utilization_percent=90.0,
    )
    update_app_settings = Mock(return_value=updated_settings)
    quantum_service = Mock()

    monkeypatch.setattr(
        settings_service_module,
        "get_current_settings",
        lambda: current_settings,
    )
    monkeypatch.setattr(settings_service_module, "update_app_settings", update_app_settings)
    monkeypatch.setattr(
        settings_service_module,
        "get_process_manager",
        lambda: SimpleNamespace(set_max_parallel_instances=Mock()),
    )
    monkeypatch.setattr(
        settings_service_module,
        "get_resource_manager",
        lambda: SimpleNamespace(update_resource_constraints=Mock()),
    )
    monkeypatch.setattr("services.get_quantum_service", lambda: quantum_service)

    service = SettingsService()

    with pytest.raises(ValidationError, match="Failed to move calculations directory"):
        service.update_settings({"calculations_directory": str(new_dir)})

    update_app_settings.assert_not_called()
    quantum_service.update_calculations_directory.assert_not_called()
    assert (old_calc_a / "metadata.json").read_text(encoding="utf-8") == "old a"
    assert (old_calc_b / "metadata.json").read_text(encoding="utf-8") == "old b"
    assert (existing_new_calc_b / "metadata.json").read_text(encoding="utf-8") == "new b"
