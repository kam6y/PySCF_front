"""Configuration manager for quantum calculation settings."""

import logging
from typing import Dict, Any, Optional, Tuple
from flask import current_app

logger = logging.getLogger(__name__)


class QuantumCalculationConfigManager:
    """
    Manager for quantum calculation configuration settings.

    This class provides convenient access to quantum calculation settings
    stored in Flask's app.config, which serves as the single source of truth
    for all application configuration.
    """

    def __init__(self):
        """
        Initialize configuration manager.

        Configuration is retrieved from Flask's app.config at runtime,
        eliminating the need for duplicate file loading and caching.
        """
        pass

    def _get_config(self) -> Dict[str, Any]:
        """
        Get quantum calculation configuration from Flask app.config.

        Returns:
            Dictionary containing quantum calculation defaults, or fallback values.
        """
        try:
            # Get configuration from Flask app.config (single source of truth)
            config = current_app.config.get('QUANTUM_CALCULATION_DEFAULTS', {})
            if config:
                return {'quantum_calculation_defaults': config}
            else:
                logger.warning("Quantum calculation defaults not found in app.config, using fallback")
                return self._get_fallback_config()
        except RuntimeError:
            # Outside Flask application context - use fallback
            logger.warning("Outside Flask context, using fallback quantum calculation configuration")
            return self._get_fallback_config()

    def _get_fallback_config(self) -> Dict[str, Any]:
        """Get fallback configuration when app.config is not available."""
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

    def reload_config(self) -> None:
        """
        Reload configuration from Flask app.config.

        Note: Since configuration is now read directly from app.config at runtime,
        this method is kept for backwards compatibility but has no effect.
        Configuration changes should be made directly to Flask's app.config.
        """
        logger.info("Configuration reload requested (configuration is now read from app.config at runtime)")
    
    def get_memory_setting(self, calculation_method: str) -> int:
        """
        Get memory setting for a specific calculation method.

        Args:
            calculation_method: The calculation method (e.g., 'CASCI', 'DFT', etc.)

        Returns:
            Memory setting in MB
        """
        config = self._get_config()
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
        config = self._get_config()
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
        config = self._get_config()
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
        config = self._get_config()
        return config["quantum_calculation_defaults"]["memory_settings"].copy()

    def get_all_cycle_settings(self) -> Dict[str, int]:
        """Get all cycle settings as a dictionary."""
        config = self._get_config()
        return config["quantum_calculation_defaults"]["cycle_settings"].copy()

    def get_all_spectrum_settings(self) -> Dict[str, Any]:
        """Get all spectrum settings as a dictionary."""
        config = self._get_config()
        return config["quantum_calculation_defaults"]["spectrum_settings"].copy()

    def get_config_info(self) -> Dict[str, Any]:
        """
        Get information about the configuration source and status.

        Returns:
            Dictionary with configuration source information.
        """
        try:
            has_config = current_app.config.get('QUANTUM_CALCULATION_DEFAULTS') is not None
            return {
                "config_source": "Flask app.config",
                "config_available": has_config,
                "using_fallback": not has_config
            }
        except RuntimeError:
            return {
                "config_source": "fallback (outside Flask context)",
                "config_available": False,
                "using_fallback": True
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