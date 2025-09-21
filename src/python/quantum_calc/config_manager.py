"""Configuration manager for quantum calculation settings."""

import json
import os
import logging
from typing import Dict, Any, Optional, List, Tuple
from pathlib import Path

logger = logging.getLogger(__name__)


class QuantumCalculationConfigManager:
    """Manager for quantum calculation configuration settings."""
    
    def __init__(self, config_path: Optional[str] = None):
        """
        Initialize configuration manager.
        
        Args:
            config_path: Path to server-config.json file. If None, auto-detects.
        """
        self.config_path = config_path or self._find_config_file()
        self._config_cache = None
        self._fallback_config = self._get_fallback_config()
        
    def _find_config_file(self) -> str:
        """Find the server-config.json file in the project hierarchy."""
        # Start from the current file's directory and walk up
        current_dir = Path(__file__).parent
        
        # Walk up the directory tree to find config/server-config.json
        for level in range(5):  # Limit search to 5 levels up
            config_path = current_dir / "config" / "server-config.json"
            if config_path.exists():
                logger.info(f"Found config file at: {config_path}")
                return str(config_path)
            current_dir = current_dir.parent
            
        # Try alternative paths
        alternative_paths = [
            os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), 'config', 'server-config.json'),
            os.path.join(os.getcwd(), 'config', 'server-config.json'),
            "config/server-config.json"
        ]
        
        for alt_path in alternative_paths:
            if os.path.exists(alt_path):
                logger.info(f"Found config file at alternative path: {alt_path}")
                return alt_path
                
        logger.warning("Could not find server-config.json, will use fallback values")
        return ""
    
    def _get_fallback_config(self) -> Dict[str, Any]:
        """Get fallback configuration when config file is not available."""
        return {
            "quantum_calculation_defaults": {
                "memory_settings": {
                    "CASCI": 6000,
                    "CASSCF": 6000,
                    "CCSD": 4000,
                    "CCSD_T": 4000,
                    "TDDFT": 4000,
                    "DFT": 2000,
                    "HF": 2000,
                    "MP2": 2000,
                    "default": 2000
                },
                "cycle_settings": {
                    "max_cycle": 150,
                    "max_cycle_macro": 50,
                    "max_cycle_micro": 4,
                    "ah_max_cycle": 30
                },
                "spectrum_settings": {
                    "ir_frequency_range": [400, 4000],
                    "frequency_threshold": 80.0
                }
            }
        }
    
    def _load_config(self) -> Dict[str, Any]:
        """Load configuration from file or return fallback."""
        if self._config_cache is not None:
            return self._config_cache
            
        if not self.config_path or not os.path.exists(self.config_path):
            logger.warning(f"Config file not found at {self.config_path}, using fallback configuration")
            self._config_cache = self._fallback_config
            return self._config_cache
            
        try:
            with open(self.config_path, 'r', encoding='utf-8') as f:
                config = json.load(f)
                
            # Merge with fallback to ensure all required keys exist
            merged_config = self._merge_with_fallback(config)
            self._config_cache = merged_config
            
            logger.info(f"Successfully loaded quantum calculation configuration from {self.config_path}")
            return self._config_cache
            
        except json.JSONDecodeError as e:
            logger.error(f"Invalid JSON in config file {self.config_path}: {e}")
            logger.warning("Using fallback configuration due to JSON parsing error")
            self._config_cache = self._fallback_config
            return self._config_cache
            
        except Exception as e:
            logger.error(f"Error loading config file {self.config_path}: {e}")
            logger.warning("Using fallback configuration due to loading error")
            self._config_cache = self._fallback_config
            return self._config_cache
    
    def _merge_with_fallback(self, config: Dict[str, Any]) -> Dict[str, Any]:
        """Merge loaded config with fallback to ensure all keys exist."""
        fallback = self._fallback_config.copy()
        
        if "quantum_calculation_defaults" in config:
            user_defaults = config["quantum_calculation_defaults"]
            fallback_defaults = fallback["quantum_calculation_defaults"]
            
            # Merge each section
            for section in ["memory_settings", "cycle_settings", "spectrum_settings"]:
                if section in user_defaults:
                    fallback_defaults[section].update(user_defaults[section])
        
        # Keep other config sections as-is
        for key, value in config.items():
            if key != "quantum_calculation_defaults":
                fallback[key] = value
                
        return fallback
    
    def reload_config(self) -> None:
        """Reload configuration from file."""
        self._config_cache = None
        self._load_config()
        logger.info("Configuration reloaded")
    
    def get_memory_setting(self, calculation_method: str) -> int:
        """
        Get memory setting for a specific calculation method.
        
        Args:
            calculation_method: The calculation method (e.g., 'CASCI', 'DFT', etc.)
            
        Returns:
            Memory setting in MB
        """
        config = self._load_config()
        memory_settings = config["quantum_calculation_defaults"]["memory_settings"]
        
        # Try exact match first
        if calculation_method in memory_settings:
            return int(memory_settings[calculation_method])
            
        # Try case-insensitive match
        method_upper = calculation_method.upper()
        if method_upper in memory_settings:
            return int(memory_settings[method_upper])
            
        # Return default
        default_memory = memory_settings.get("default", 2000)
        logger.debug(f"No specific memory setting found for {calculation_method}, using default: {default_memory} MB")
        return int(default_memory)
    
    def get_cycle_setting(self, setting_name: str) -> int:
        """
        Get cycle setting value.
        
        Args:
            setting_name: Name of the cycle setting (e.g., 'max_cycle', 'max_cycle_macro')
            
        Returns:
            Cycle setting value
        """
        config = self._load_config()
        cycle_settings = config["quantum_calculation_defaults"]["cycle_settings"]
        
        if setting_name in cycle_settings:
            return int(cycle_settings[setting_name])
            
        # Fallback values for unknown settings
        fallback_values = {
            "max_cycle": 150,
            "max_cycle_macro": 50,
            "max_cycle_micro": 4,
            "ah_max_cycle": 30
        }
        
        default_value = fallback_values.get(setting_name, 150)
        logger.debug(f"No setting found for {setting_name}, using default: {default_value}")
        return default_value
    
    def get_spectrum_setting(self, setting_name: str) -> Any:
        """
        Get spectrum setting value.
        
        Args:
            setting_name: Name of the spectrum setting
            
        Returns:
            Spectrum setting value
        """
        config = self._load_config()
        spectrum_settings = config["quantum_calculation_defaults"]["spectrum_settings"]
        
        if setting_name in spectrum_settings:
            return spectrum_settings[setting_name]
            
        # Fallback values
        fallback_values = {
            "ir_frequency_range": [400, 4000],
            "frequency_threshold": 80.0
        }
        
        default_value = fallback_values.get(setting_name)
        logger.debug(f"No spectrum setting found for {setting_name}, using default: {default_value}")
        return default_value
    
    def get_ir_frequency_range(self) -> Tuple[float, float]:
        """Get IR frequency range as a tuple."""
        range_list = self.get_spectrum_setting("ir_frequency_range")
        if isinstance(range_list, list) and len(range_list) >= 2:
            return (float(range_list[0]), float(range_list[1]))
        return (400.0, 4000.0)
    
    def get_frequency_threshold(self) -> float:
        """Get frequency threshold for vibrational analysis."""
        return float(self.get_spectrum_setting("frequency_threshold"))
    
    def get_all_memory_settings(self) -> Dict[str, int]:
        """Get all memory settings as a dictionary."""
        config = self._load_config()
        return config["quantum_calculation_defaults"]["memory_settings"].copy()
    
    def get_all_cycle_settings(self) -> Dict[str, int]:
        """Get all cycle settings as a dictionary."""
        config = self._load_config()
        return config["quantum_calculation_defaults"]["cycle_settings"].copy()
    
    def get_all_spectrum_settings(self) -> Dict[str, Any]:
        """Get all spectrum settings as a dictionary."""
        config = self._load_config()
        return config["quantum_calculation_defaults"]["spectrum_settings"].copy()
    
    def get_config_info(self) -> Dict[str, Any]:
        """Get information about the configuration source and status."""
        return {
            "config_path": self.config_path,
            "config_exists": os.path.exists(self.config_path) if self.config_path else False,
            "using_fallback": self._config_cache == self._fallback_config if self._config_cache else True,
            "cache_loaded": self._config_cache is not None
        }


# Global instance for easy access
_config_manager = None

def get_config_manager() -> QuantumCalculationConfigManager:
    """Get the global configuration manager instance."""
    global _config_manager
    if _config_manager is None:
        _config_manager = QuantumCalculationConfigManager()
    return _config_manager

def reload_global_config() -> None:
    """Reload the global configuration."""
    global _config_manager
    if _config_manager is not None:
        _config_manager.reload_config()

# Convenience functions for common operations
def get_memory_for_method(calculation_method: str) -> int:
    """Convenience function to get memory setting for a calculation method."""
    return get_config_manager().get_memory_setting(calculation_method)

def get_max_cycle() -> int:
    """Convenience function to get max_cycle setting."""
    return get_config_manager().get_cycle_setting("max_cycle")

def get_max_cycle_macro() -> int:
    """Convenience function to get max_cycle_macro setting."""
    return get_config_manager().get_cycle_setting("max_cycle_macro")

def get_max_cycle_micro() -> int:
    """Convenience function to get max_cycle_micro setting."""
    return get_config_manager().get_cycle_setting("max_cycle_micro")

def get_ah_max_cycle() -> int:
    """Convenience function to get ah_max_cycle setting."""
    return get_config_manager().get_cycle_setting("ah_max_cycle")

def get_ir_range() -> Tuple[float, float]:
    """Convenience function to get IR frequency range."""
    return get_config_manager().get_ir_frequency_range()

def get_spectrum_setting(setting_name: str) -> Any:
    """Convenience function to get spectrum setting."""
    return get_config_manager().get_spectrum_setting(setting_name)

def get_frequency_threshold() -> float:
    """Convenience function to get frequency threshold."""
    return get_config_manager().get_frequency_threshold()