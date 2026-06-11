"""
Settings management service.

This service encapsulates application settings management logic,
providing a unified interface for both API endpoints and AI agent tools.
"""

import logging
from typing import Dict, Any

from quantum_calc import (
    get_process_manager,
    get_current_settings,
    update_app_settings,
    CalculationDirectoryMigration,
    mask_settings,
)
from quantum_calc.resource_manager import get_resource_manager
from .exceptions import ServiceError, ValidationError

logger = logging.getLogger(__name__)


class SettingsService:
    """Service for application settings management.

    Note: ``get_settings`` and ``update_settings`` intentionally return the
    **real** settings values (including ``gemini_api_key``) because internal
    callers (e.g. the Gemini agent chat in ``api/agent.py``) need the
    plaintext key.  Masking for HTTP responses is handled at the API
    boundary (``api/settings.py::_mask_api_key_for_response``).

    Security: the API key is currently stored as plaintext JSON at rest.
    See the threat-model note in ``quantum_calc.settings_manager`` for the
    rationale and the deferred keyring-based fix.
    """

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
    def _rollback_successful_directory_migration(
        move_result: Dict[str, Any],
        current_calc_dir: str,
    ) -> None:
        """Move calculation data back after a later settings update failure."""
        if (
            move_result.get("success") is not True
            or move_result.get("moved_count", 0) <= 0
        ):
            return

        new_path = move_result.get("new_path")
        old_path = move_result.get("old_path", current_calc_dir)
        if not new_path or not old_path:
            logger.warning(
                "Cannot roll back calculations directory migration because paths are missing: "
                f"{move_result}"
            )
            return

        logger.warning(
            "Settings save failed after calculations directory migration; "
            f"rolling back {new_path} -> {old_path}"
        )
        try:
            rollback_migration = CalculationDirectoryMigration(base_dir=new_path)
            rollback_result = rollback_migration.move_calculations_directory(old_path)
        except Exception as rollback_error:
            logger.error(
                f"Failed to roll back calculations directory migration: {rollback_error}",
                exc_info=True,
            )
            raise ServiceError(
                "Failed to update settings after moving calculations directory, "
                f"and rollback failed: {rollback_error}"
            ) from rollback_error

        if rollback_result.get("success") is not True:
            logger.error(
                f"Failed to roll back calculations directory migration: {rollback_result}"
            )
            raise ServiceError(
                "Failed to update settings after moving calculations directory, "
                f"and rollback failed: {rollback_result.get('message', 'Move operation failed')}"
            )

        logger.info(
            "Rolled back calculations directory migration after settings save failure"
        )

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

            logger.info("Successfully retrieved settings")
            return settings.model_dump(mode="json")

        except Exception as e:
            logger.error(f"Failed to retrieve settings: {e}", exc_info=True)
            raise ServiceError(f"Failed to retrieve settings: {str(e)}")

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

            # Preserve stored values for fields that are easily lost on
            # round-trip saves.  The GET endpoint masks gemini_api_key to "",
            # so a naive save-back would wipe it; research_email can likewise
            # arrive as None when the caller omits it (the Pydantic model
            # defaults Optional fields to None, and model_dump always emits
            # the key).  Only a non-empty string is treated as an intentional
            # update; absent / None / "" values are stripped so the stored
            # value is preserved.  Build a new dict to honour the immutability
            # coding principle.
            _PRESERVE_WHEN_EMPTY = ("gemini_api_key", "research_email")
            keys_to_strip = {k for k in _PRESERVE_WHEN_EMPTY if not new_settings.get(k)}
            if keys_to_strip:
                new_settings = {
                    k: v for k, v in new_settings.items() if k not in keys_to_strip
                }

            # Get current settings to detect changes
            current_settings = get_current_settings()

            # Check if calculations_directory is changing
            new_calc_dir = new_settings.get("calculations_directory")
            current_calc_dir = current_settings.calculations_directory

            move_result = None
            if new_calc_dir and new_calc_dir != current_calc_dir:
                logger.info(
                    f"Calculations directory changing from {current_calc_dir} to {new_calc_dir}"
                )
                self._raise_if_calculations_directory_change_blocked()

                try:
                    migration = CalculationDirectoryMigration(base_dir=current_calc_dir)

                    # Move calculations to new directory
                    move_result = migration.move_calculations_directory(new_calc_dir)

                    if move_result.get("success") is not True:
                        logger.warning(
                            f"Calculations directory migration failed: {move_result}"
                        )
                        self._raise_for_failed_migration(move_result)

                    logger.info(
                        f"Successfully moved calculations: {move_result['message']}"
                    )

                except ValueError as move_error:
                    logger.error(
                        f"Invalid calculations directory migration: {move_error}"
                    )
                    raise ValidationError(
                        f"Failed to move calculations directory: {str(move_error)}"
                    )
                except OSError as move_error:
                    logger.error(
                        f"Failed to move calculations directory: {move_error}",
                        exc_info=True,
                    )
                    raise ServiceError(
                        f"Failed to move calculations directory: {str(move_error)}"
                    )
                except ServiceError:
                    raise
                except Exception as move_error:
                    logger.error(
                        f"Failed to move calculations directory: {move_error}",
                        exc_info=True,
                    )
                    raise ServiceError(
                        f"Failed to move calculations directory: {str(move_error)}"
                    )

            # Update settings
            try:
                updated_settings = update_app_settings(new_settings)
            except Exception:
                if new_calc_dir and new_calc_dir != current_calc_dir and move_result:
                    self._rollback_successful_directory_migration(
                        move_result,
                        current_calc_dir,
                    )
                raise

            # Update process manager with new parallel instance limit
            try:
                process_manager = get_process_manager()
                process_manager.set_max_parallel_instances(
                    updated_settings.max_parallel_instances
                )
            except Exception as pm_error:
                logger.warning(f"Failed to update process manager settings: {pm_error}")

            # Update resource manager with new resource constraints
            try:
                resource_manager = get_resource_manager()
                resource_manager.update_resource_constraints(
                    max_cpu_utilization_percent=updated_settings.max_cpu_utilization_percent,
                    max_memory_utilization_percent=updated_settings.max_memory_utilization_percent,
                )
            except Exception as rm_error:
                logger.warning(
                    f"Failed to update resource manager settings: {rm_error}"
                )

            # Update quantum service with new calculations directory if it changed
            if new_calc_dir and new_calc_dir != current_calc_dir:
                try:
                    # Import here to avoid circular dependency
                    from . import get_quantum_service

                    quantum_service = get_quantum_service()
                    quantum_service.update_calculations_directory(
                        updated_settings.calculations_directory
                    )
                    logger.info(
                        "Updated QuantumService with new calculations directory"
                    )
                except Exception as qs_error:
                    logger.warning(
                        f"Failed to update quantum service settings: {qs_error}"
                    )

            result = updated_settings.model_dump(mode="json")

            logger.info("Successfully updated settings")
            return result

        except ValueError as e:
            logger.error(f"Invalid settings values: {e}")
            raise ValidationError(f"Invalid settings: {str(e)}")
        except ServiceError:
            raise
        except Exception as e:
            logger.error(f"Failed to update settings: {e}", exc_info=True)
            raise ServiceError(f"Failed to update settings: {str(e)}")
