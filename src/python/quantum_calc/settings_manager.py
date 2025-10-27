"""Application settings manager for PySCF Native App."""

import os
import json
import multiprocessing
import logging
from typing import Dict, Any, Optional
from pathlib import Path
from generated_models import AppSettings

logger = logging.getLogger(__name__)

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
        
        self.settings_file.parent.mkdir(parents=True, exist_ok=True)
        
        logger.info(f"Settings file: {self.settings_file}")
    
    def get_default_settings(self) -> AppSettings:
        """Get default application settings."""
        try:
            import psutil
            total_cores = psutil.cpu_count(logical=True)
            total_memory_mb = int(psutil.virtual_memory().total / (1024 * 1024))
        except ImportError:
            # Fallback when psutil is not available
            total_cores = multiprocessing.cpu_count()
            total_memory_mb = 4096  # Conservative 4GB estimate

        # Default calculations directory (with PySCF_calculations subfolder)
        default_calc_dir = str(Path.home() / "PySCF_calculations")

        return AppSettings(
            max_parallel_instances=min(4, total_cores),
            max_cpu_utilization_percent=95.0,
            max_memory_utilization_percent=95.0,
            system_total_cores=total_cores,
            system_total_memory_mb=total_memory_mb,
            calculations_directory=default_calc_dir,
            gemini_api_key=None,
            tavily_api_key=None,
            research_email="pyscf-research-agent@example.com"
        )
    
    def load_settings(self) -> AppSettings:
        """
        Load settings from file with automatic migration support.
        
        Returns:
            AppSettings: Loaded settings, with migration applied if necessary.
        """
        try:
            if self.settings_file.exists():
                with open(self.settings_file, 'r', encoding='utf-8') as f:
                    data = json.load(f)
                
                try:
                    # Try to validate loaded data using Pydantic model
                    settings = AppSettings(**data)
                    logger.info(f"Loaded settings: {settings}")
                    return settings
                except Exception as validation_error:
                    # Migration needed - merge existing data with defaults
                    logger.warning(f"Settings validation failed: {validation_error}. Performing migration.")
                    return self._migrate_settings(data)
            else:
                # Create new settings file with defaults
                default_settings = self.get_default_settings()
                self.save_settings(default_settings)
                logger.info(f"Created new settings file with defaults: {default_settings}")
                return default_settings
                
        except (json.JSONDecodeError, ValueError, TypeError) as e:
            logger.warning(f"Failed to parse settings from {self.settings_file}: {e}. Resetting to defaults.")
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
                        logger.warning(f"Settings migration: ignoring invalid type for {key}: {type(value)} (expected {expected_type})")
            
            # Create validated settings object
            migrated_settings = AppSettings(**merged_data)
            
            # Save migrated settings back to file
            if self.save_settings(migrated_settings):
                logger.info(f"Successfully migrated settings: {migrated_settings}")
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
            settings_dict = settings.model_dump()
            
            # Write to temporary file first, then rename for atomic operation
            temp_file = self.settings_file.with_suffix('.json.tmp')
            
            with open(temp_file, 'w', encoding='utf-8') as f:
                json.dump(settings_dict, f, indent=2, ensure_ascii=False)
            
            # Atomic rename
            temp_file.replace(self.settings_file)
            
            logger.info(f"Saved settings: {settings}")
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