"""Application settings manager for PySCF_front.

Threat model – credential storage at rest
==========================================
The Gemini API key (gemini_api_key) is persisted as plaintext inside the
JSON settings file (~/.pyscf_native_app/app-settings.json).  Current
mitigations:

* File permissions hardened to 0600 (owner-only) on POSIX systems.
* Directory permissions hardened to 0700.
* HTTP responses never expose the raw key (masked at the API boundary).
* Log output masks the key via mask_settings().

Residual risk: any process running under the same OS user account can read
the settings file.  For a single-user local desktop application this is an
accepted tradeoff — the same user already has access to process memory,
browser cookies and other per-user secrets.

Full fix (deferred): store the key in the OS credential store (macOS
Keychain / Windows Credential Manager / freedesktop Secret Service) via
the keyring library.  This requires adding keyring to the packaged
conda environment and updating save_settings / load_settings to read/write
the key through the keyring API instead of the JSON file.
"""

import os
import json
import multiprocessing
import logging
from typing import Dict, Any, Optional
from pathlib import Path
from generated_models import AppSettings, Timezone

logger = logging.getLogger(__name__)

_SENSITIVE_KEYS = {"gemini_api_key", "research_email"}

# Restrictive file-system permission modes for sensitive data
_DIR_MODE = 0o700  # Owner-only: rwx------
_FILE_MODE = 0o600  # Owner-only: rw-------


def _secure_permissions(path: Path, mode: int) -> None:
    """Set restrictive file-system permissions on a path (POSIX only).

    Silently skips on non-POSIX platforms (e.g. Windows) where
    ``os.chmod`` has different semantics.

    Args:
        path: File or directory path to secure.
        mode: Octal permission mode (e.g. 0o700 for directories).
    """
    if os.name != "posix":
        return
    try:
        os.chmod(path, mode)
    except OSError as exc:
        logger.warning(
            "Failed to set permissions %s on %s: %s",
            oct(mode),
            path,
            exc,
        )


def mask_settings(settings) -> dict:
    """Return settings as a dict with sensitive fields masked for logging."""
    if hasattr(settings, "model_dump"):
        d = settings.model_dump(mode="json")
    elif isinstance(settings, dict):
        d = dict(settings)
    else:
        return {"<unserializable>": str(type(settings))}
    for key in _SENSITIVE_KEYS:
        if key in d and d[key] is not None:
            d[key] = "***"
    return d


class SettingsManager:
    """Manager for application settings persistence and retrieval."""

    def __init__(self, settings_file: Optional[str] = None):
        """
        Initialize the settings manager.

        Args:
            settings_file: Path to the settings file. If None, uses default location.
        """
        if settings_file is None:
            # Default settings file location in the user data directory
            base_dir = Path(os.path.expanduser("~"))
            app_data_dir = base_dir / ".pyscf_native_app"
            app_data_dir.mkdir(exist_ok=True)
            self.settings_file = app_data_dir / "app-settings.json"
        else:
            self.settings_file = Path(settings_file)

        # F3: Create directory with restrictive mode from the start.
        # os.makedirs mode= only sets the leaf directory's perms (parents
        # use umask), which is acceptable — only the leaf holds sensitive
        # files.  The post-create chmod below is an idempotent safety net
        # for directories that already existed with lax permissions.
        os.makedirs(str(self.settings_file.parent), mode=_DIR_MODE, exist_ok=True)

        # Idempotent: harden directory permissions for existing installs
        _secure_permissions(self.settings_file.parent, _DIR_MODE)

        # F3: Harden existing settings file permissions if it already exists
        if self.settings_file.exists():
            _secure_permissions(self.settings_file, _FILE_MODE)

        logger.info(f"Settings file: {self.settings_file}")

    def get_default_settings(self) -> AppSettings:
        """Get default application settings."""

        def fallback_cpu_count() -> int:
            try:
                return max(1, int(multiprocessing.cpu_count() or 1))
            except (NotImplementedError, ValueError, TypeError):
                return 1

        try:
            import psutil

            total_cores = psutil.cpu_count(logical=True) or fallback_cpu_count()
            total_memory_mb = int(psutil.virtual_memory().total / (1024 * 1024))
        except ImportError:
            # Fallback when psutil is not available
            total_cores = fallback_cpu_count()
            total_memory_mb = 4096  # Conservative 4GB estimate
        except Exception as exc:
            logger.warning(
                f"Failed to detect system resources via psutil: {exc}. Using conservative defaults."
            )
            total_cores = fallback_cpu_count()
            total_memory_mb = 4096

        total_cores = max(1, int(total_cores or 1))
        total_memory_mb = max(1, int(total_memory_mb or 4096))

        # Default calculations directory (with PySCF_calculations subfolder)
        default_calc_dir = str(Path.home() / "PySCF_calculations")

        return AppSettings(
            max_parallel_instances=min(4, total_cores),
            max_cpu_utilization_percent=95.0,
            max_memory_utilization_percent=95.0,
            system_total_cores=total_cores,
            system_total_memory_mb=total_memory_mb,
            calculations_directory=default_calc_dir,
            timezone=Timezone.UTC,
            gemini_api_key=None,
            research_email="pyscf-research-agent@example.com",
            gpu_acceleration_enabled=False,
        )

    def load_settings(self) -> AppSettings:
        """
        Load settings from file with automatic migration support.

        Returns:
            AppSettings: Loaded settings, with migration applied if necessary.
        """
        try:
            if self.settings_file.exists():
                with open(self.settings_file, "r", encoding="utf-8") as f:
                    data = json.load(f)

                try:
                    # Try to validate loaded data using Pydantic model
                    settings = AppSettings(**data)
                    logger.info(f"Loaded settings: {mask_settings(settings)}")
                    return settings
                except Exception as validation_error:
                    # Migration needed - merge existing data with defaults
                    logger.warning(
                        f"Settings validation failed: {validation_error}. Performing migration."
                    )
                    return self._migrate_settings(data)
            else:
                # Create new settings file with defaults
                default_settings = self.get_default_settings()
                self.save_settings(default_settings)
                logger.info(
                    f"Created new settings file with defaults: {mask_settings(default_settings)}"
                )
                return default_settings

        except (json.JSONDecodeError, ValueError, TypeError) as e:
            logger.warning(
                f"Failed to parse settings from {self.settings_file}: {e}. Resetting to defaults."
            )
            default_settings = self.get_default_settings()
            self.save_settings(default_settings)
            return default_settings
        except Exception as e:
            logger.error(f"Unexpected error loading settings: {e}. Using defaults.")
            return self.get_default_settings()

    def _migrate_settings(self, existing_data: Dict[str, Any]) -> AppSettings:
        """
        Migrate existing settings data to current schema.

        Args:
            existing_data: Existing settings data from file.

        Returns:
            AppSettings: Migrated settings.
        """
        try:
            # Get default settings as baseline
            default_settings = self.get_default_settings()
            default_dict = default_settings.model_dump()

            # Merge existing data with defaults (existing data takes priority where valid)
            merged_data = default_dict.copy()

            # Update with existing valid values
            for key, value in existing_data.items():
                if key in default_dict:
                    # Validate the type matches expected type
                    expected_type = type(default_dict[key])
                    if isinstance(value, expected_type):
                        merged_data[key] = value
                    else:
                        logger.warning(
                            f"Settings migration: ignoring invalid type for {key}: {type(value)} (expected {expected_type})"
                        )

            # Create validated settings object
            migrated_settings = AppSettings(**merged_data)

            # Save migrated settings back to file
            if self.save_settings(migrated_settings):
                logger.info(
                    f"Successfully migrated settings: {mask_settings(migrated_settings)}"
                )
            else:
                logger.warning("Failed to save migrated settings to file")

            return migrated_settings

        except Exception as e:
            logger.error(f"Settings migration failed: {e}. Using defaults.")
            return self.get_default_settings()

    def save_settings(self, settings: AppSettings) -> bool:
        """
        Save settings to file.

        Args:
            settings: Settings to save.

        Returns:
            bool: True if save was successful, False otherwise.
        """
        try:
            # Convert Pydantic model to dict for JSON serialization
            # mode='json' ensures enums are serialized as their values
            settings_dict = settings.model_dump(mode="json")

            # Write to temporary file first, then rename for atomic operation
            temp_file = self.settings_file.with_suffix(".json.tmp")

            # F1: On POSIX, create the temp file with restrictive permissions
            # from the start via os.open() so there is no world-readable
            # window.  On non-POSIX (Windows), fall back to plain open().
            if os.name == "posix":
                fd = os.open(
                    str(temp_file),
                    os.O_CREAT | os.O_WRONLY | os.O_TRUNC,
                    _FILE_MODE,
                )
                with os.fdopen(fd, "w", encoding="utf-8") as f:
                    json.dump(settings_dict, f, indent=2, ensure_ascii=False)
            else:
                with open(temp_file, "w", encoding="utf-8") as f:
                    json.dump(settings_dict, f, indent=2, ensure_ascii=False)

            # Idempotent safety net: ensure permissions are correct even if
            # the temp file previously existed with lax permissions.
            _secure_permissions(temp_file, _FILE_MODE)

            # Atomic rename (preserves inode permissions from the temp file)
            temp_file.replace(self.settings_file)

            logger.info(f"Saved settings: {mask_settings(settings)}")
            return True

        except Exception as e:
            logger.error(f"Failed to save settings to {self.settings_file}: {e}")
            return False

    def update_settings(self, updates: Dict[str, Any]) -> AppSettings:
        """
        Update specific settings fields.

        Args:
            updates: Dictionary of field updates.

        Returns:
            AppSettings: Updated settings.

        Raises:
            ValueError: If invalid settings are provided.
        """
        current_settings = self.load_settings()
        current_dict = current_settings.model_dump()

        # Apply updates
        current_dict.update(updates)

        # Validate updated settings
        updated_settings = AppSettings(**current_dict)

        # Save updated settings
        if self.save_settings(updated_settings):
            return updated_settings
        else:
            raise RuntimeError("Failed to save updated settings")


# Global settings manager instance
_settings_manager: Optional[SettingsManager] = None


def get_settings_manager() -> SettingsManager:
    """Get the global settings manager instance."""
    global _settings_manager
    if _settings_manager is None:
        _settings_manager = SettingsManager()
    return _settings_manager


def get_current_settings() -> AppSettings:
    """Get current application settings."""
    return get_settings_manager().load_settings()


def update_app_settings(updates: Dict[str, Any]) -> AppSettings:
    """
    Update application settings.

    Args:
        updates: Dictionary of settings updates.

    Returns:
        Updated settings.
    """
    return get_settings_manager().update_settings(updates)
