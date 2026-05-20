"""
Unit tests for SettingsManager defaults.
"""

import pytest
import quantum_calc.settings_manager as settings_manager
from quantum_calc.settings_manager import SettingsManager


def test_default_settings_gpu_disabled(tmp_path):
    """
    GIVEN SettingsManager default settings
    WHEN get_default_settings is called
    THEN GPU acceleration is disabled by default
    """
    manager = SettingsManager(settings_file=str(tmp_path / "settings.json"))
    settings = manager.get_default_settings()

    assert settings.gpu_acceleration_enabled is False


def test_default_settings_handles_unknown_logical_cpu_count(tmp_path, monkeypatch):
    """
    GIVEN psutil cannot determine logical CPU count
    WHEN default settings are generated
    THEN settings still contain a usable positive core count
    """
    psutil = pytest.importorskip("psutil")
    monkeypatch.setattr(psutil, "cpu_count", lambda logical=True: None)
    monkeypatch.setattr(settings_manager.multiprocessing, "cpu_count", lambda: 2)

    manager = SettingsManager(settings_file=str(tmp_path / "settings.json"))
    settings = manager.get_default_settings()

    assert settings.system_total_cores == 2
    assert settings.max_parallel_instances == 2
