"""
Unit tests for SettingsManager defaults and file-system security.
"""

import os
import stat

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


# ---------------------------------------------------------------------------
# File-system permission tests (POSIX only)
# ---------------------------------------------------------------------------



@pytest.mark.skipif(os.name != "posix", reason="Permission tests require POSIX")
class TestSettingsFilePermissions:
    """Verify that SettingsManager creates files/dirs with restrictive modes."""

    def test_save_settings_file_mode_is_0600(self, tmp_path):
        """
        GIVEN a SettingsManager pointed at a temp directory
        WHEN settings are saved
        THEN the settings file has mode 0o600 (owner read/write only)
        """
        # Arrange
        settings_file = tmp_path / "app_data" / "settings.json"
        manager = SettingsManager(settings_file=str(settings_file))

        # Act
        default_settings = manager.get_default_settings()
        manager.save_settings(default_settings)

        # Assert
        actual_mode = stat.S_IMODE(os.stat(settings_file).st_mode)
        assert actual_mode == 0o600, f"Settings file mode {oct(actual_mode)} != 0o600"

    def test_settings_directory_mode_is_0700(self, tmp_path):
        """
        GIVEN a SettingsManager pointed at a temp directory
        WHEN the manager is initialised (directory auto-created)
        THEN the settings directory has mode 0o700 (owner rwx only)
        """
        # Arrange
        settings_dir = tmp_path / "secure_dir"
        settings_file = settings_dir / "settings.json"

        # Act
        SettingsManager(settings_file=str(settings_file))

        # Assert
        actual_mode = stat.S_IMODE(os.stat(settings_dir).st_mode)
        assert (
            actual_mode == 0o700
        ), f"Settings directory mode {oct(actual_mode)} != 0o700"

    def test_permissions_explicit_not_inherited_from_umask(self, tmp_path):
        """
        GIVEN a permissive umask (0o000)
        WHEN SettingsManager saves settings
        THEN the file is still 0o600 and the directory is 0o700
              (proving permissions are set explicitly, not inherited)
        """
        # Arrange
        settings_dir = tmp_path / "umask_test"
        settings_file = settings_dir / "settings.json"
        old_umask = os.umask(0)

        try:
            # Act
            manager = SettingsManager(settings_file=str(settings_file))
            manager.save_settings(manager.get_default_settings())

            # Assert — directory
            dir_mode = stat.S_IMODE(os.stat(settings_dir).st_mode)
            assert (
                dir_mode == 0o700
            ), f"Directory mode {oct(dir_mode)} != 0o700 under umask(0)"

            # Assert — file
            file_mode = stat.S_IMODE(os.stat(settings_file).st_mode)
            assert (
                file_mode == 0o600
            ), f"File mode {oct(file_mode)} != 0o600 under umask(0)"
        finally:
            os.umask(old_umask)
