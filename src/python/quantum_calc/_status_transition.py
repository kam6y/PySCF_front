"""Calculation status transitions and notification management."""

import logging
from enum import Enum
from typing import Optional, Callable

logger = logging.getLogger(__name__)


class CalculationStatus(Enum):
    """Calculation lifecycle states."""
    WAITING = "waiting"
    RUNNING = "running"
    PAUSING = "pausing"
    PAUSED = "paused"
    COMPLETED = "completed"
    ERROR = "error"


class CalculationStatusManager:
    """Manages calculation status transitions, file persistence, and realtime update notifications."""

    def __init__(self, notification_callback: Optional[Callable] = None):
        self.notification_callback = notification_callback

    def transition(self, calculation_id: str, new_status: CalculationStatus,
                   error_message: Optional[str] = None) -> None:
        """
        Execute a status transition: persist to file + send realtime update notification.

        This is used by the parent process (CalculationProcessManager) for all status
        updates. Worker processes write status directly via CalculationRepository since they
        cannot access parent process objects.
        """
        from quantum_calc._calculation_repository import CalculationRepository
        from quantum_calc import get_current_settings

        settings = get_current_settings()
        repository = CalculationRepository(base_dir=settings.calculations_directory)
        calc_dir = str(repository.resolve_calculation_path(calculation_id))

        status_str = new_status.value
        repository.save_calculation_status(calc_dir, status_str)

        if new_status == CalculationStatus.ERROR and error_message:
            repository.save_calculation_results(calc_dir, {'error': error_message})
            logger.info(f"Calculation {calculation_id} transitioned to 'error': {error_message}")
        else:
            logger.info(f"Calculation {calculation_id} transitioned to '{status_str}'")

        self.notify(calculation_id, status_str, error_message)

    def notify(self, calculation_id: str, status: str,
               error_message: Optional[str] = None) -> None:
        """Send realtime update notification for a calculation status change."""
        if self.notification_callback is None:
            logger.debug(f"Realtime update notification not available for calculation {calculation_id}")
            return

        try:
            self.notification_callback(calculation_id, status, error_message)
            logger.debug(f"Sent realtime update notification for calculation {calculation_id} with status {status}")
        except Exception as e:
            logger.warning(f"Failed to send realtime update notification for calculation {calculation_id}: {e}")
