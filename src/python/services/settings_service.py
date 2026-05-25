"""
Settings management service.

This service encapsulates application settings management logic,
providing a unified interface for both API endpoints and AI agent tools.
"""

import logging
from typing import Dict, Any

from quantum_calc import (
    get_process_manager, get_current_settings, update_app_settings,
    CalculationDirectoryMigration, mask_settings,
)
from quantum_calc.resource_manager import get_resource_manager
from .exceptions import ServiceError, ValidationError

logger = logging.getLogger(__name__)


class SettingsService:
    """Service for application settings management."""

    DIRECTORY_CHANGE_BLOCKED_MESSAGE = (
        "計算中またはキュー中は計算ディレクトリを変更できません。"
        "計算が完了してから再度お試しください。"
    )

    @staticmethod
    def _raise_for_failed_migration(move_result: Dict[str, Any]) -> None:
        """Convert a failed directory migration result into a service exception."""
        message = move_result.get("message", "Move operation failed")
        failed_moves = move_result.get("failed_moves", [])
        rollback_errors = move_result.get("rollback_errors", [])

        is_user_fixable_conflict = (
            failed_moves
            and not rollback_errors
            and all(
                failed_move.get("reason") == "Destination already exists"
                for failed_move in failed_moves
            )
        )

        error_message = f"Failed to move calculations directory: {message}"
        if is_user_fixable_conflict:
            raise ValidationError(error_message)
        raise ServiceError(error_message)

    @staticmethod
    def _raise_if_calculations_directory_change_blocked() -> None:
        """Reject directory changes while calculations are managed in memory."""
        process_manager = get_process_manager()
        active_calculations = process_manager.get_active_calculations() or []
        queued_calculations = process_manager.get_queued_calculations() or []

        if active_calculations or queued_calculations:
            raise ValidationError(SettingsService.DIRECTORY_CHANGE_BLOCKED_MESSAGE)
    
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
            return settings.model_dump(mode='json')
            
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
            logger.info(f"Updating application settings: {mask_settings(new_settings)}")

            # Get current settings to detect changes
            current_settings = get_current_settings()

            # Check if calculations_directory is changing
            new_calc_dir = new_settings.get('calculations_directory')
            current_calc_dir = current_settings.calculations_directory

            move_result = None
            if new_calc_dir and new_calc_dir != current_calc_dir:
                logger.info(f"Calculations directory changing from {current_calc_dir} to {new_calc_dir}")
                self._raise_if_calculations_directory_change_blocked()

                try:
                    migration = CalculationDirectoryMigration(base_dir=current_calc_dir)

                    # Move calculations to new directory
                    move_result = migration.move_calculations_directory(new_calc_dir)

                    if move_result.get('success') is not True:
                        logger.warning(f"Calculations directory migration failed: {move_result}")
                        self._raise_for_failed_migration(move_result)

                    logger.info(f"Successfully moved calculations: {move_result['message']}")

                except ValueError as move_error:
                    logger.error(f"Invalid calculations directory migration: {move_error}")
                    raise ValidationError(
                        f'Failed to move calculations directory: {str(move_error)}'
                    )
                except OSError as move_error:
                    logger.error(f"Failed to move calculations directory: {move_error}", exc_info=True)
                    raise ServiceError(f'Failed to move calculations directory: {str(move_error)}')
                except ServiceError:
                    raise
                except Exception as move_error:
                    logger.error(f"Failed to move calculations directory: {move_error}", exc_info=True)
                    raise ServiceError(f'Failed to move calculations directory: {str(move_error)}')

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

            # Update quantum service with new calculations directory if it changed
            if new_calc_dir and new_calc_dir != current_calc_dir:
                try:
                    # Import here to avoid circular dependency
                    from . import get_quantum_service
                    quantum_service = get_quantum_service()
                    quantum_service.update_calculations_directory(updated_settings.calculations_directory)
                    logger.info(f"Updated QuantumService with new calculations directory")
                except Exception as qs_error:
                    logger.warning(f"Failed to update quantum service settings: {qs_error}")

            result = updated_settings.model_dump(mode='json')

            logger.info(f"Successfully updated settings")
            return result

        except ValueError as e:
            logger.error(f"Invalid settings values: {e}")
            raise ValidationError(f'Invalid settings: {str(e)}')
        except ServiceError:
            raise
        except Exception as e:
            logger.error(f"Failed to update settings: {e}", exc_info=True)
            raise ServiceError(f'Failed to update settings: {str(e)}')
