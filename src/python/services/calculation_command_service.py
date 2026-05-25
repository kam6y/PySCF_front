"""State-changing calculation service operations."""

import logging
import os
import shutil
from contextlib import suppress
from typing import Any

from quantum_calc import (
    FileManagerError,
    GeometryError,
    InputError,
    ProcessManagerError,
    get_current_settings,
)

from .calculation_service_context import CalculationServiceContext
from .exceptions import (
    InsufficientResourcesError,
    NotFoundError,
    PermissionDeniedError,
    ResourceUnavailableError,
    ServiceError,
    ValidationError,
)

logger = logging.getLogger(__name__)


class CalculationCommandService:
    """Handles calculation submission and lifecycle mutations."""

    GPU4PYSCF_SUPPORTED_METHODS = frozenset({"DFT", "HF", "TDDFT"})

    def __init__(
        self,
        context: CalculationServiceContext,
        query_service: Any,
    ) -> None:
        self.context = context
        self.query_service = query_service

    def _with_runtime_settings_snapshot(self, params: dict[str, Any]) -> dict[str, Any]:
        """Return calculation parameters with runtime settings fixed for this job."""
        settings = get_current_settings()
        return {
            **params,
            "gpu_acceleration_enabled": bool(
                getattr(settings, "gpu_acceleration_enabled", False)
            ),
        }

    def validate_calculation_parameters(self, params: dict[str, Any]) -> str | None:
        """
        Validate calculation parameters for compatibility and theoretical correctness.

        Returns:
            None if validation passes, error message string if validation fails.
        """
        calculation_method = params.get("calculation_method")

        from quantum_calc.method_defaults import validate_parameters_for_method

        is_valid, applicability_error = validate_parameters_for_method(
            calculation_method,
            params,
        )
        if not is_valid:
            return applicability_error

        gpu_acceleration_enabled = params.get("gpu_acceleration_enabled")
        if gpu_acceleration_enabled is None:
            settings = get_current_settings()
            gpu_acceleration_enabled = getattr(
                settings,
                "gpu_acceleration_enabled",
                False,
            )
        if (
            bool(gpu_acceleration_enabled)
            and calculation_method not in self.GPU4PYSCF_SUPPORTED_METHODS
        ):
            supported_methods = ", ".join(sorted(self.GPU4PYSCF_SUPPORTED_METHODS))
            return (
                "GPU acceleration is enabled, but "
                f"{calculation_method} is not supported by GPU4PySCF. "
                f"Supported GPU methods: {supported_methods}. "
                "Disable GPU acceleration or choose a supported method."
            )

        if calculation_method in {"DFT", "TDDFT"} and not params.get(
            "exchange_correlation"
        ):
            return (
                f"{calculation_method} method requires an exchange-correlation "
                "functional to be specified"
            )

        if calculation_method == "TDDFT" and (
            not params.get("tddft_nstates") or params.get("tddft_nstates") < 1
        ):
            return (
                "TDDFT method requires tddft_nstates to be specified and greater "
                "than 0"
            )

        if calculation_method in ["CASCI", "CASSCF"]:
            if not params.get("ncas") or params.get("ncas") < 1:
                return (
                    f"{calculation_method} method requires ncas (active space "
                    "orbitals) to be specified and greater than 0"
                )
            if not params.get("nelecas") or params.get("nelecas") < 1:
                return (
                    f"{calculation_method} method requires nelecas (active space "
                    "electrons) to be specified and greater than 0"
                )
            if params.get("nelecas") > 2 * params.get("ncas"):
                return (
                    f"{calculation_method} method: nelecas ({params.get('nelecas')}) "
                    f"cannot exceed 2 * ncas ({2 * params.get('ncas')})"
                )

        if params.get("spin", 0) < 0:
            return "Spin multiplicity (2S) cannot be negative"

        if params.get("charges") is not None and abs(params.get("charges")) > 10:
            logger.warning(
                f"High molecular charge ({params.get('charges')}) - please verify "
                "this is correct"
            )

        return None

    def start_calculation(self, params: dict[str, Any]) -> dict[str, Any]:
        """
        Start a new quantum chemistry calculation.

        Raises:
            ValidationError: If parameters are invalid.
            ResourceUnavailableError: If process manager is unavailable.
            InsufficientResourcesError: If system resources are insufficient.
            ServiceError: For other failures.
        """
        try:
            params = self._with_runtime_settings_snapshot(params)

            validation_error = self.validate_calculation_parameters(params)
            if validation_error:
                logger.warning(f"Parameter validation failed: {validation_error}")
                raise ValidationError(f"Invalid parameters: {validation_error}")

            try:
                process_manager = self.context.get_process_manager()
                self.context.recover_stale_non_terminal_calculations(process_manager)
            except ProcessManagerError as e:
                logger.error(f"Process manager error: {e}")
                raise ResourceUnavailableError(
                    "System initialization error: Unable to initialize calculation "
                    "system. Please check system resources and try again."
                )
            except Exception as submit_error:
                error_message = f"Failed to submit calculation: {str(submit_error)}"
                logger.error(
                    f"Unexpected error during calculation submission: {submit_error}"
                )
                raise ServiceError(error_message)

            try:
                calc_dir = self.context.repository.create_calculation_dir(params["name"])
                calculation_id = os.path.basename(calc_dir)

                self.context.repository.save_calculation_parameters(calc_dir, params)
                logger.info(
                    "Created calculation directory and saved parameters for "
                    f"calculation {calculation_id}"
                )
            except Exception as file_error:
                logger.error(f"Failed to set up calculation files: {file_error}")
                raise ServiceError(f"Failed to initialize calculation: {str(file_error)}")

            try:
                try:
                    from quantum_calc import update_process_manager_settings

                    update_process_manager_settings()
                except Exception as settings_error:
                    logger.warning(
                        f"Failed to update process manager settings: {settings_error}"
                    )

                logger.info(f"About to submit calculation {calculation_id}")
                success, initial_status, waiting_reason = (
                    process_manager.submit_calculation(calculation_id, params)
                )
                logger.info(
                    "Submit result: "
                    f"success={success}, status={initial_status}, "
                    f"reason={waiting_reason}"
                )

                if not success:
                    error_message = (
                        waiting_reason
                        if waiting_reason
                        else "Failed to submit calculation to process pool."
                    )
                    self.context.repository.save_calculation_status(calc_dir, "error")
                    self.context.repository.save_calculation_results(
                        calc_dir,
                        {"error": error_message},
                    )
                    logger.error(
                        f"Failed to submit calculation {calculation_id} to process "
                        f"pool: {error_message}"
                    )

                    error_instance = self.context.build_calculation_instance(
                        calculation_id,
                        params,
                        "error",
                    )
                    error_instance["error"] = error_message
                    return error_instance

                current_status, current_waiting_reason = (
                    self.context.repository.read_calculation_status_details(calc_dir)
                )
                if current_status in self.context.TERMINAL_STATUSES:
                    logger.info(
                        "Calculation %s already reached terminal status %s during submit",
                        calculation_id,
                        current_status,
                    )
                    return self.context.build_calculation_instance(
                        calculation_id,
                        params,
                        current_status,
                        current_waiting_reason,
                    )

                self.context.repository.save_calculation_status(
                    calc_dir,
                    initial_status,
                    waiting_reason,
                )
                logger.info(
                    f"Queued calculation {calculation_id} for molecule '{params['name']}'"
                )

                return self.context.build_calculation_instance(
                    calculation_id,
                    params,
                    initial_status,
                    waiting_reason,
                )

            except ProcessManagerError as e:
                with suppress(Exception):
                    shutil.rmtree(calc_dir, ignore_errors=True)
                logger.error(f"Process manager error: {e}")
                raise ResourceUnavailableError(
                    "System initialization error: Unable to initialize calculation "
                    "system. Please check system resources and try again."
                )
            except Exception as submit_error:
                self.context.repository.save_calculation_status(calc_dir, "error")
                error_message = f"Failed to submit calculation: {str(submit_error)}"
                self.context.repository.save_calculation_results(
                    calc_dir,
                    {"error": error_message},
                )
                logger.error(
                    f"Unexpected error during calculation submission: {submit_error}"
                )
                raise ServiceError(error_message)

        except (InputError, GeometryError) as e:
            logger.warning(f"Invalid calculation parameters: {e}")
            raise ValidationError(f"Invalid input parameters: {str(e)}")
        except FileManagerError as e:
            logger.error(f"File management error during calculation setup: {e}")
            raise ServiceError("Failed to set up calculation files.")
        except OSError as e:
            logger.error(f"System error during calculation setup: {e}")
            raise InsufficientResourcesError(
                "Insufficient system resources to start calculation."
            )
        except PermissionError as e:
            logger.error(f"Permission error during calculation setup: {e}")
            raise PermissionDeniedError(
                "System permission error. Please contact administrator."
            )

    def update_calculation(
        self,
        calculation_id: str,
        new_name: str,
    ) -> dict[str, Any]:
        """
        Update calculation metadata.

        Raises:
            NotFoundError: If calculation not found.
            ValidationError: If inputs are invalid.
            ServiceError: For file access failures.
        """
        try:
            self.context.resolve_calculation_path(calculation_id)
            result_id = self.context.repository.rename_calculation(
                calculation_id,
                new_name,
            )
            if not result_id:
                raise NotFoundError(f'Calculation "{calculation_id}" not found.')

            logger.info(
                f"Updated display name for calculation {calculation_id} to '{new_name}'"
            )
            return {
                "message": "Calculation renamed successfully.",
                "name": new_name,
            }
        except FileManagerError as e:
            logger.error(f"File manager error updating calculation {calculation_id}: {e}")
            raise ServiceError("Unable to update calculation data.")
        except ValueError as e:
            logger.warning(f"Invalid input for calculation {calculation_id}: {e}")
            raise ValidationError("Invalid input data.")
        except OSError as e:
            logger.error(f"System error updating calculation {calculation_id}: {e}")
            raise ServiceError("System error updating calculation.")

    def delete_calculation(self, calculation_id: str) -> dict[str, Any]:
        """
        Delete a calculation and its files.

        Raises:
            NotFoundError: If calculation not found.
            ValidationError: If calculation cannot be deleted in its current state.
            ServiceError: For file access failures.
        """
        try:
            calc_path = self.context.resolve_calculation_path(calculation_id)
            self.context.recover_stale_non_terminal_calculations()

            if not os.path.isdir(calc_path):
                raise NotFoundError(f'Calculation "{calculation_id}" not found.')

            status, _ = self.context.repository.read_calculation_status_details(
                calc_path
            )
            non_deletable_statuses = {"pending", "running", "waiting", "pausing"}
            if status in non_deletable_statuses:
                logger.warning(
                    "Cannot delete calculation %s with non-terminal status %s",
                    calculation_id,
                    status,
                )
                raise ValidationError(
                    f'Cannot delete calculation "{calculation_id}" while it is '
                    f"{status}. Please pause or wait for the calculation to complete "
                    "first."
                )

            process_manager = self.context.get_process_manager()

            if process_manager.is_running(calculation_id):
                logger.warning(f"Cannot delete running calculation {calculation_id}")
                raise ValidationError(
                    f'Cannot delete calculation "{calculation_id}" while it is '
                    "running. Please pause or wait for the calculation to complete "
                    "first."
                )

            shutil.rmtree(calc_path)
            logger.info(f"Deleted calculation directory: {calc_path}")

            return {
                "message": (
                    f'Calculation "{calculation_id}" has been deleted successfully'
                ),
                "deleted_id": calculation_id,
            }

        except ProcessManagerError as e:
            logger.error(f"Process manager error deleting calculation {calculation_id}: {e}")
            raise ResourceUnavailableError("Process manager unavailable.")
        except FileManagerError as e:
            logger.error(f"File manager error deleting calculation {calculation_id}: {e}")
            raise ServiceError("Unable to access calculation files.")
        except PermissionError as e:
            logger.error(f"Permission error deleting calculation {calculation_id}: {e}")
            raise PermissionDeniedError("Permission denied deleting calculation.")
        except OSError as e:
            logger.error(f"System error deleting calculation {calculation_id}: {e}")
            raise ServiceError("System error deleting calculation files.")
        except ValueError as e:
            logger.warning(f"Invalid calculation ID for deletion {calculation_id}: {e}")
            raise ValidationError("Invalid calculation ID.")

    def pause_calculation(self, calculation_id: str) -> dict[str, Any]:
        """
        Pause a running calculation.

        Raises:
            NotFoundError: If calculation not found.
            ValidationError: If calculation is not pausable.
            ServiceError: If pause fails.
        """
        try:
            logger.info(f"Pausing calculation: {calculation_id}")
            calc_path = self.context.resolve_calculation_path(calculation_id)
            if not os.path.isdir(calc_path):
                raise NotFoundError(f'Calculation "{calculation_id}" not found.')

            process_manager = self.context.get_process_manager()
            self.context.recover_stale_non_terminal_calculations(process_manager)

            status, _ = self.context.repository.read_calculation_status_details(
                calc_path
            )
            if status != "running":
                raise ValidationError(f"Calculation is not running (status: {status})")

            success = process_manager.pause_calculation(calculation_id)

            if not success:
                raise ServiceError("Failed to pause calculation")

            logger.info(f"Calculation pause requested successfully: {calculation_id}")

            return {
                "message": (
                    "Pause request accepted. Calculation will pause after current "
                    "iteration."
                ),
                "calculation_id": calculation_id,
            }

        except ValueError as e:
            message = str(e)
            logger.error(f"Cannot pause calculation {calculation_id}: {message}")
            if "not found" in message.lower():
                raise NotFoundError(f'Calculation "{calculation_id}" not found.') from e
            raise ValidationError(message) from e
        except (NotFoundError, ValidationError):
            raise
        except Exception as e:
            logger.error(f"Error pausing calculation {calculation_id}: {e}", exc_info=True)
            raise ServiceError(f"Failed to pause calculation: {str(e)}")

    def resume_calculation(self, calculation_id: str) -> dict[str, Any]:
        """
        Resume a paused calculation.

        Raises:
            NotFoundError: If calculation not found.
            ValidationError: If calculation cannot be resumed.
            ServiceError: If resume fails.
        """
        try:
            logger.info(f"Resuming calculation: {calculation_id}")
            calc_path = self.context.resolve_calculation_path(calculation_id)
            if not os.path.isdir(calc_path):
                raise NotFoundError(f'Calculation "{calculation_id}" not found.')

            process_manager = self.context.get_process_manager()
            process_manager.resume_calculation(calculation_id)

            logger.info(f"Calculation resumed successfully: {calculation_id}")

            calculation_details = self.query_service.get_calculation_details(
                calculation_id
            )

            return {
                "message": "Calculation resumed from checkpoint",
                "calculation": calculation_details["calculation"],
            }

        except ValueError as e:
            message = str(e)
            logger.error(f"Cannot resume calculation {calculation_id}: {message}")
            if "not found" in message.lower():
                raise NotFoundError(f'Calculation "{calculation_id}" not found.') from e
            raise ValidationError(message) from e
        except (NotFoundError, ValidationError):
            raise
        except Exception as e:
            logger.error(
                f"Error resuming calculation {calculation_id}: {e}",
                exc_info=True,
            )
            raise ServiceError(f"Failed to resume calculation: {str(e)}")
