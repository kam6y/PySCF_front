"""
Pause Manager for quantum chemistry calculations.

This module provides functionality to pause and resume calculations gracefully.
"""

import threading
from pathlib import Path
from typing import Set
import logging

logger = logging.getLogger(__name__)


class PauseManager:
    """
    Manages pause/resume requests for calculations.

    This class provides a thread-safe mechanism to:
    - Request pauses for running calculations
    - Check if a pause is requested
    - Create/check flag files for inter-process communication
    """

    def __init__(self):
        """Initialize the pause manager."""
        self._pause_requests: Set[str] = set()
        self._lock = threading.Lock()

    def request_pause(self, calc_id: str) -> None:
        """
        Request pause for a calculation.

        Args:
            calc_id: Unique calculation ID
        """
        with self._lock:
            self._pause_requests.add(calc_id)
            logger.info(f"Pause requested for calculation: {calc_id}")

    def is_pause_requested(self, calc_id: str) -> bool:
        """
        Check if pause is requested for a calculation.

        Args:
            calc_id: Unique calculation ID

        Returns:
            True if pause is requested, False otherwise
        """
        with self._lock:
            return calc_id in self._pause_requests

    def clear_pause_request(self, calc_id: str) -> None:
        """
        Clear pause request after handling.

        Args:
            calc_id: Unique calculation ID
        """
        with self._lock:
            self._pause_requests.discard(calc_id)
            logger.info(f"Pause request cleared for calculation: {calc_id}")

    def create_pause_flag_file(self, calc_dir: str) -> None:
        """
        Create a flag file to signal pause to worker process.

        This file is used for inter-process communication when the worker
        is running in a separate process.

        Args:
            calc_dir: Calculation directory path
        """
        flag_file = Path(calc_dir) / ".pause_requested"
        try:
            flag_file.touch()
            logger.debug(f"Created pause flag file: {flag_file}")
        except Exception as e:
            logger.error(f"Failed to create pause flag file: {e}")

    def check_pause_flag_file(self, calc_dir: str) -> bool:
        """
        Check if pause flag file exists.

        Args:
            calc_dir: Calculation directory path

        Returns:
            True if flag file exists, False otherwise
        """
        flag_file = Path(calc_dir) / ".pause_requested"
        return flag_file.exists()

    def remove_pause_flag_file(self, calc_dir: str) -> None:
        """
        Remove pause flag file.

        Args:
            calc_dir: Calculation directory path
        """
        flag_file = Path(calc_dir) / ".pause_requested"
        try:
            if flag_file.exists():
                flag_file.unlink()
                logger.debug(f"Removed pause flag file: {flag_file}")
        except Exception as e:
            logger.error(f"Failed to remove pause flag file: {e}")


# Global instance
pause_manager = PauseManager()
