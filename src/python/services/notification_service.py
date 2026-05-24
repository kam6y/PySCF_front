"""WebSocket notification service."""
import logging
from typing import Optional
from datetime import datetime

from websocket.event_loop_bridge import schedule_coroutine

logger = logging.getLogger(__name__)


class NotificationService:
    """Service for sending WebSocket notifications to connected clients."""

    def __init__(self):
        self._socketio = None
        logger.info("NotificationService initialized (socketio not yet bound)")

    def bind_socketio(self, socketio) -> None:
        """Bind SocketIO instance to this service."""
        self._socketio = socketio
        logger.info("Socket.IO AsyncServer bound to NotificationService")

    def send_calculation_update(
        self,
        calculation_id: str,
        status: str,
        error_message: Optional[str] = None
    ) -> None:
        """Send immediate WebSocket notification for calculation status changes."""
        if self._socketio is None:
            logger.warning(
                "Cannot send notification for %s: Socket.IO not bound",
                calculation_id,
            )
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
                logger.warning("Calculation directory not found: %s", calc_dir)
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

        except Exception:
            logger.exception("Error building WebSocket notification payload")
            return

        scheduled_future = schedule_coroutine(
            self._socketio.emit(
                "calculation_update",
                calculation_instance,
                room="global_updates",
            )
        )
        if scheduled_future is None:
            logger.warning(
                "Dropped WebSocket notification for %s because scheduling failed",
                calculation_id,
            )
            return

        logger.debug(
            "Scheduled WebSocket notification for %s with status %s",
            calculation_id,
            status,
        )


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
    logger.info("Global NotificationService bound to Socket.IO")
