"""
Settings management service.

This service encapsulates application settings management logic,
providing a unified interface for both API endpoints and AI agent tools.
"""

import logging
from typing import Dict, Any

from quantum_calc import get_process_manager, get_current_settings, update_app_settings
from quantum_calc.resource_manager import get_resource_manager
from .exceptions import ServiceError, ValidationError

logger = logging.getLogger(__name__)


class SettingsService:
    """Service for application settings management."""
    
    def get_settings(self) -> Dict[str, Any]:
        """
        Get current application settings.
        
        Returns:
            Dict containing current settings
            
        Raises:
            ServiceError: If retrieval fails
        """
        try:
            logger.info("Getting application settings")
            
            settings = get_current_settings()
            
            logger.info(f"Successfully retrieved settings")
            return settings.model_dump()
            
        except Exception as e:
            logger.error(f"Failed to retrieve settings: {e}", exc_info=True)
            raise ServiceError(f'Failed to retrieve settings: {str(e)}')
    
    def update_settings(self, new_settings: Dict[str, Any]) -> Dict[str, Any]:
        """
        Update application settings.
        
        Args:
            new_settings: Dictionary of new settings values
            
        Returns:
            Dict containing updated settings
            
        Raises:
            ValidationError: If settings values are invalid
            ServiceError: For other errors
        """
        try:
            logger.info(f"Updating application settings: {new_settings}")
            
            # Update settings
            updated_settings = update_app_settings(new_settings)
            
            # Update process manager with new parallel instance limit
            try:
                process_manager = get_process_manager()
                process_manager.set_max_parallel_instances(updated_settings.max_parallel_instances)
            except Exception as pm_error:
                logger.warning(f"Failed to update process manager settings: {pm_error}")
            
            # Update resource manager with new resource constraints
            try:
                resource_manager = get_resource_manager()
                resource_manager.update_resource_constraints(
                    max_cpu_utilization_percent=updated_settings.max_cpu_utilization_percent,
                    max_memory_utilization_percent=updated_settings.max_memory_utilization_percent
                )
            except Exception as rm_error:
                logger.warning(f"Failed to update resource manager settings: {rm_error}")
            
            logger.info(f"Successfully updated settings")
            return updated_settings.model_dump()
            
        except ValueError as e:
            logger.error(f"Invalid settings values: {e}")
            raise ValidationError(f'Invalid settings: {str(e)}')
        except Exception as e:
            logger.error(f"Failed to update settings: {e}", exc_info=True)
            raise ServiceError(f'Failed to update settings: {str(e)}')
