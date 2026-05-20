from pathlib import Path
from types import SimpleNamespace
from typing import Any
from unittest.mock import Mock

import pytest

import services.settings_service as settings_service_module
from services.exceptions import ValidationError
from services.settings_service import SettingsService


DIRECTORY_CHANGE_BLOCKED_MESSAGE = "計算中またはキュー中は計算ディレクトリを変更できません"


class _UpdatedSettings(SimpleNamespace):
    def model_dump(self, mode: str | None = None) -> dict[str, Any]:
        return self.__dict__


def _settings_update_dependencies(
    old_dir: Path,
    new_dir: Path,
    active_calculations: list[str] | None = None,
    queued_calculations: list[str] | None = None,
) -> tuple[SimpleNamespace, _UpdatedSettings, Mock, Mock, Mock, Mock]:
    current_settings = SimpleNamespace(calculations_directory=str(old_dir))
    updated_settings = _UpdatedSettings(
        calculations_directory=str(new_dir),
        max_parallel_instances=2,
        max_cpu_utilization_percent=90.0,
        max_memory_utilization_percent=90.0,
    )
    process_manager = Mock()
    process_manager.get_active_calculations.return_value = active_calculations or []
    process_manager.get_queued_calculations.return_value = queued_calculations or []
    update_app_settings = Mock(return_value=updated_settings)
    quantum_service = Mock()
    migration_class = Mock()
    migration = migration_class.return_value
    migration.move_calculations_directory.return_value = {
        "success": True,
        "message": "Moved calculations directory",
    }
    return (
        current_settings,
        updated_settings,
        process_manager,
        update_app_settings,
        quantum_service,
        migration_class,
    )


@pytest.mark.parametrize(
    ("active_calculations", "queued_calculations"),
    [
        (["active-calc"], []),
        ([], ["queued-calc"]),
    ],
)
def test_update_settings_rejects_directory_migration_with_active_or_queued_calculations(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    active_calculations: list[str],
    queued_calculations: list[str],
) -> None:
    """
    GIVEN active or queued calculations exist
    WHEN settings update requests a calculations_directory change
    THEN directory migration and settings updates are rejected before side effects
    """
    old_dir = tmp_path / "old_calculations"
    new_dir = tmp_path / "new_calculations"
    (
        current_settings,
        _updated_settings,
        process_manager,
        update_app_settings,
        quantum_service,
        migration_class,
    ) = _settings_update_dependencies(
        old_dir,
        new_dir,
        active_calculations=active_calculations,
        queued_calculations=queued_calculations,
    )

    monkeypatch.setattr(
        settings_service_module,
        "get_current_settings",
        lambda: current_settings,
    )
    monkeypatch.setattr(settings_service_module, "get_process_manager", lambda: process_manager)
    monkeypatch.setattr(settings_service_module, "update_app_settings", update_app_settings)
    monkeypatch.setattr(
        settings_service_module,
        "CalculationDirectoryMigration",
        migration_class,
    )
    monkeypatch.setattr("services.get_quantum_service", lambda: quantum_service)

    service = SettingsService()

    with pytest.raises(ValidationError, match=DIRECTORY_CHANGE_BLOCKED_MESSAGE):
        service.update_settings({"calculations_directory": str(new_dir)})

    migration_class.assert_not_called()
    migration_class.return_value.move_calculations_directory.assert_not_called()
    update_app_settings.assert_not_called()
    quantum_service.update_calculations_directory.assert_not_called()


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
        lambda: SimpleNamespace(
            get_active_calculations=Mock(return_value=[]),
            get_queued_calculations=Mock(return_value=[]),
            set_max_parallel_instances=Mock(),
        ),
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
