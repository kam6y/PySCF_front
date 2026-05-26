"""Calculation update notification service."""

import logging
import os
from typing import Optional

from websocket.event_loop_bridge import schedule_coroutine

logger = logging.getLogger(__name__)


class NotificationService:
    """Service for publishing calculation update notifications."""

    def send_calculation_update(
        self,
        calculation_id: str,
        status: str,
        error_message: Optional[str] = None,
    ) -> None:
        """Publish immediate notification for calculation status changes."""
        try:
            from quantum_calc import CalculationRepository, get_current_settings
            from services.calculation_update_payload import build_calculation_instance
            from services.calculation_update_stream import (
                get_calculation_update_stream_hub,
            )

            settings = get_current_settings()
            file_manager = CalculationRepository(
                base_dir=settings.calculations_directory
            )
            calc_dir = os.path.join(file_manager.get_base_directory(), calculation_id)

            if not os.path.exists(calc_dir):
                logger.warning("Calculation directory not found: %s", calc_dir)
                return

            calculation_instance = build_calculation_instance(
                calculation_id,
                calc_dir,
                file_manager,
                status_override=status,
                error_message=error_message,
            )

        except Exception:
            logger.exception("Error building calculation notification payload")
            return

        scheduled_future = schedule_coroutine(
            get_calculation_update_stream_hub().publish_calculation_update(
                calculation_instance
            )
        )
        if scheduled_future is None:
            logger.warning(
                "Dropped calculation notification for %s because scheduling failed",
                calculation_id,
            )
            return

        logger.debug(
            "Scheduled calculation notification for %s with status %s",
            calculation_id,
            status,
        )


_notification_service: Optional[NotificationService] = None


def get_notification_service() -> NotificationService:
    global _notification_service
    if _notification_service is None:
        _notification_service = NotificationService()
    return _notification_service
