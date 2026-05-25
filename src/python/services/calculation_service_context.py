"""Shared context and helpers for calculation service boundaries."""

import logging
from dataclasses import dataclass
from typing import Any

from quantum_calc import (
    CalculationRepository,
    CubeArtifactService,
    get_current_settings,
    get_process_manager,
)

from .exceptions import ValidationError

logger = logging.getLogger(__name__)


@dataclass
class CalculationServiceContext:
    """Shared calculation storage and cross-boundary helpers."""

    repository: CalculationRepository
    cube_service: CubeArtifactService

    TERMINAL_STATUSES = frozenset({"completed", "error", "paused"})
    NON_TERMINAL_STATUSES = frozenset({"pending", "running", "waiting", "pausing"})
    RESTART_INTERRUPTED_MESSAGE = (
        "Calculation interrupted because the backend was restarted."
    )

    @classmethod
    def from_settings(cls) -> "CalculationServiceContext":
        """Build a shared context from the current application settings."""
        settings = get_current_settings()
        calculations_dir = settings.calculations_directory
        return cls(
            repository=CalculationRepository(base_dir=calculations_dir),
            cube_service=CubeArtifactService(base_dir=calculations_dir),
        )

    def set_base_directory(self, new_directory: str) -> None:
        """Update all calculation storage helpers to use the same base directory."""
        self.repository.set_base_directory(new_directory)
        self.cube_service.set_base_directory(new_directory)

    def get_process_manager(self) -> Any:
        """Return the shared calculation process manager."""
        return get_process_manager()

    def build_calculation_instance(
        self,
        calculation_id: str,
        parameters: dict[str, Any],
        status: str,
        waiting_reason: str | None = None,
    ) -> dict[str, Any]:
        """Build the calculation instance shape returned by public APIs."""
        instance = {
            "id": calculation_id,
            "name": parameters["name"],
            "status": status,
            "createdAt": parameters["created_at"],
            "updatedAt": parameters["created_at"],
            "parameters": parameters,
        }

        if waiting_reason is not None:
            instance["waitingReason"] = waiting_reason

        return instance

    def resolve_calculation_path(self, calculation_id: str) -> str:
        """Resolve a calculation ID from an external request into a safe path."""
        try:
            return str(self.repository.resolve_calculation_path(calculation_id))
        except ValueError as e:
            raise ValidationError(str(e)) from e

    def recover_stale_non_terminal_calculations(
        self,
        process_manager: Any | None = None,
    ) -> None:
        """Mark persisted non-terminal calculations as error when no worker owns them."""
        try:
            manager = process_manager or self.get_process_manager()
            active_ids = set(manager.get_active_calculations() or [])
            queued_ids = set(manager.get_queued_calculations() or [])
        except Exception as e:
            logger.warning(f"Skipping stale calculation recovery: {e}")
            return

        managed_ids = active_ids | queued_ids
        for calculation in self.repository.list_calculations():
            calculation_id = calculation.get("id")
            if not calculation_id or calculation_id in managed_ids:
                continue

            try:
                calc_dir = str(self.repository.resolve_calculation_path(calculation_id))
                status, _ = self.repository.read_calculation_status_details(calc_dir)
                if status not in self.NON_TERMINAL_STATUSES:
                    continue

                logger.warning(
                    "Recovering stale calculation %s from status %s to error",
                    calculation_id,
                    status,
                )
                existing_results = self.repository.read_calculation_results(calc_dir)
                self.repository.save_calculation_status(calc_dir, "error")
                if existing_results is None:
                    self.repository.save_calculation_results(
                        calc_dir,
                        {"error": self.RESTART_INTERRUPTED_MESSAGE},
                    )
            except Exception as e:
                logger.warning(
                    "Failed to recover stale calculation %s: %s",
                    calculation_id,
                    e,
                )
