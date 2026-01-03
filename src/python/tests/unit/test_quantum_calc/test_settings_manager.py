"""
Unit tests for SettingsManager defaults.
"""

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
