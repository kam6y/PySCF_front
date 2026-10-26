"""Read-only calculation service operations."""

import logging
import os
from datetime import datetime
from typing import Any

from quantum_calc import (
    FileManagerError,
    ProcessManagerError,
    get_all_supported_parameters,
)

from .calculation_service_context import CalculationServiceContext
from .exceptions import (
    NotFoundError,
    PermissionDeniedError,
    ResourceUnavailableError,
    ServiceError,
)

logger = logging.getLogger(__name__)


class CalculationQueryService:
    """Handles supported-parameter, listing, detail, and system-status queries."""

    def __init__(self, context: CalculationServiceContext) -> None:
        self.context = context

    def get_supported_parameters(self) -> dict[str, Any]:
        """
        Get supported quantum chemistry parameters.

        Raises:
            ServiceError: If parameter retrieval fails.
        """
        try:
            logger.info("Getting supported quantum chemistry parameters")
            parameters = get_all_supported_parameters()
            logger.info("Successfully retrieved supported parameters")
            return parameters
        except Exception as e:
            logger.error(f"Error getting supported parameters: {e}", exc_info=True)
            raise ServiceError(f"Failed to retrieve supported parameters: {str(e)}")

    def list_calculations(
        self,
        name_query: str | None = None,
        status: str | None = None,
        calculation_method: str | None = None,
        basis_function: str | None = None,
        date_from: str | None = None,
        date_to: str | None = None,
    ) -> dict[str, Any]:
        """
        List available calculations with optional filtering.

        Raises:
            ServiceError: If listing fails.
        """
        try:
            self.context.recover_stale_non_terminal_calculations()

            calculations = self.context.repository.list_calculations(
                name_query=name_query,
                status=status,
                calculation_method=calculation_method,
                basis_function=basis_function,
                date_from=date_from,
                date_to=date_to,
            )

            return {
                "calculations": calculations,
                "count": len(calculations),
            }
        except FileManagerError as e:
            logger.error(f"File manager error while listing calculations: {e}")
            raise ServiceError("Unable to access calculation directory.")
        except PermissionError as e:
            logger.error(f"Permission error while listing calculations: {e}")
            raise PermissionDeniedError(
                "Permission denied accessing calculation directory."
            )
        except OSError as e:
            logger.error(f"System error while listing calculations: {e}")
            raise ServiceError("System error accessing calculation files.")

    def get_calculation_details(self, calculation_id: str) -> dict[str, Any]:
        """
        Get detailed information about a specific calculation.

        Raises:
            NotFoundError: If calculation not found.
            ServiceError: For file access failures.
        """
        try:
            calc_path = self.context.resolve_calculation_path(calculation_id)
            self.context.recover_stale_non_terminal_calculations()

            if not os.path.isdir(calc_path):
                raise NotFoundError(f'Calculation "{calculation_id}" not found.')

            parameters = (
                self.context.repository.read_calculation_parameters(calc_path) or {}
            )
            results = self.context.repository.read_calculation_results(calc_path)
            status, waiting_reason = (
                self.context.repository.read_calculation_status_details(calc_path)
            )

            display_name = self.context.repository.get_display_name(parameters)
            creation_date = parameters.get(
                "created_at",
                datetime.fromtimestamp(os.path.getmtime(calc_path)).isoformat(),
            )

            calculation_instance = {
                "id": calculation_id,
                "name": display_name,
                "status": status,
                "createdAt": creation_date,
                "updatedAt": datetime.fromtimestamp(
                    os.path.getmtime(calc_path)
                ).isoformat(),
                "workingDirectory": calc_path,
                "parameters": parameters,
                "results": results,
            }

            if waiting_reason is not None:
                calculation_instance["waitingReason"] = waiting_reason

            return {
                "calculation": calculation_instance,
                "files": {
                    "checkpoint_exists": self.context.repository.file_exists(
                        calc_path,
                        "calculation.chk",
                    ),
                    "parameters_file_exists": parameters is not None,
                    "results_file_exists": results is not None,
                },
            }
        except FileManagerError as e:
            logger.error(f"File manager error getting calculation {calculation_id}: {e}")
            raise ServiceError("Unable to access calculation data.")
        except OSError as e:
            logger.error(f"System error getting calculation {calculation_id}: {e}")
            raise ServiceError("System error accessing calculation files.")

    def get_calculation_status(self) -> dict[str, Any]:
        """
        Get status information about the calculation system.

        Raises:
            ResourceUnavailableError: If process manager is unavailable.
            ServiceError: If the process manager is not properly initialized.
        """
        try:
            process_manager = self.context.get_process_manager()
            active_calculations = process_manager.get_active_calculations()

            return {
                "process_pool": {
                    "max_workers": process_manager.max_workers,
                    "active_calculations": active_calculations,
                    "active_count": len(active_calculations),
                    "is_shutdown": process_manager._shutdown,
                },
                "system": {
                    "cpu_count": os.cpu_count(),
                },
            }
        except ProcessManagerError as e:
            logger.error(f"Process manager error getting status: {e}")
            raise ResourceUnavailableError("Process manager is unavailable.")
        except AttributeError as e:
            logger.error(f"Process manager attribute error: {e}")
            raise ServiceError("Process manager not properly initialized.")
