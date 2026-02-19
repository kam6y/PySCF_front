"""WebSocket notification service."""
import logging
from typing import Optional
from datetime import datetime

logger = logging.getLogger(__name__)


class NotificationService:
    """Service for sending WebSocket notifications to connected clients."""

    def __init__(self):
        self._socketio = None
        logger.info("NotificationService initialized (socketio not yet bound)")

    def bind_socketio(self, socketio) -> None:
        """Bind SocketIO instance to this service."""
        self._socketio = socketio
        logger.info("SocketIO instance bound to NotificationService")

    def send_calculation_update(
        self,
        calculation_id: str,
        status: str,
        error_message: Optional[str] = None
    ) -> None:
        """Send immediate WebSocket notification for calculation status changes."""
        if self._socketio is None:
            logger.warning(f"Cannot send notification for {calculation_id}: SocketIO not bound")
            return

        try:
            from quantum_calc import get_current_settings, CalculationRepository
            import os

            settings = get_current_settings()
            file_manager = CalculationRepository(
                base_dir=settings.calculations_directory
            )
            calc_dir = os.path.join(
                file_manager.get_base_directory(),
                calculation_id
            )

            if not os.path.exists(calc_dir):
                logger.warning(f"Calculation directory not found: {calc_dir}")
                return

            # Read current data
            parameters = file_manager.read_calculation_parameters(calc_dir) or {}
            results = file_manager.read_calculation_results(calc_dir)
            display_name = file_manager.get_display_name(calculation_id, parameters)

            # Build calculation instance
            calculation_instance = {
                'id': calculation_id,
                'name': display_name,
                'status': status,
                'createdAt': parameters.get('created_at', datetime.now().isoformat()),
                'updatedAt': datetime.now().isoformat(),
                'parameters': parameters,
                'results': results,
                'workingDirectory': calc_dir,
            }

            if error_message:
                calculation_instance['error'] = error_message
                calculation_instance['errorMessage'] = error_message

            # Send to global updates room
            # (calculation-specific room is handled by file_watcher via handlers.py)
            self._socketio.emit(
                'calculation_update',
                calculation_instance,
                room='global_updates'
            )

            logger.debug(f"Sent WebSocket notification for {calculation_id} with status {status}")

        except Exception as e:
            logger.error(f"Error sending WebSocket notification: {e}", exc_info=True)


# Global singleton
_notification_service: Optional[NotificationService] = None


def get_notification_service() -> NotificationService:
    """Get the global NotificationService instance."""
    global _notification_service
    if _notification_service is None:
        _notification_service = NotificationService()
    return _notification_service


def bind_notification_service(socketio) -> None:
    """Bind SocketIO instance to the global NotificationService."""
    service = get_notification_service()
    service.bind_socketio(socketio)
    logger.info("Global NotificationService bound to SocketIO")
