"""Calculation queue management and resource-based scheduling."""

import logging
import threading
from datetime import datetime
from typing import Optional, Callable, List, Dict
from dataclasses import dataclass
from concurrent.futures import Future

from .config_manager import get_memory_for_method
from .resource_manager import AllocationStatus

logger = logging.getLogger(__name__)


@dataclass
class QueuedCalculation:
    """Represents a calculation waiting in the queue."""
    calculation_id: str
    parameters: dict
    created_at: datetime
    waiting_reason: Optional[str] = None

    def __lt__(self, other):
        """For priority queue ordering by creation time."""
        return self.created_at < other.created_at


class CalculationQueueScheduler:
    """Manages the calculation queue and resource-based scheduling decisions."""

    def __init__(self, resource_manager, max_parallel_instances: int):
        self.resource_manager = resource_manager
        self.max_parallel_instances = max_parallel_instances
        self.calculation_queue: List[QueuedCalculation] = []
        self._queue_lock = threading.Lock()
        self._queue_processing = False

    def check_and_decide(self, parameters: dict, active_count: int,
                         max_slots: int) -> tuple:
        """
        Check resource availability and decide whether to start, queue, or reject.
        Returns (AllocationStatus, reason).
        """
        user_cpu_cores = parameters.get('cpu_cores') or 1
        user_memory_mb = parameters.get('memory_mb') or get_memory_for_method(
            parameters.get('calculation_method', 'DFT'))
        calculation_method = parameters.get('calculation_method', 'DFT')

        if self.resource_manager is None:
            return AllocationStatus.CAN_START, "Resource checking disabled"

        try:
            return self.resource_manager.check_allocation_status(
                cpu_cores=user_cpu_cores,
                memory_mb=user_memory_mb,
                calculation_method=calculation_method,
                active_count=active_count,
                max_slots=max_slots
            )
        except Exception as e:
            logger.warning(f"Resource allocation check failed: {e}. Proceeding with CAN_START.")
            return AllocationStatus.CAN_START, "Resource checking failed - proceeding without constraints"

    def enqueue(self, calculation_id: str, parameters: dict, reason: str) -> None:
        """Add a calculation to the waiting queue."""
        queued_calc = QueuedCalculation(
            calculation_id=calculation_id,
            parameters=parameters,
            created_at=datetime.now(),
            waiting_reason=reason
        )
        self.calculation_queue.append(queued_calc)
        self.calculation_queue.sort(key=lambda x: x.created_at)
        logger.info(f"Added calculation {calculation_id} to queue (position {len(self.calculation_queue)}): {reason}")

    def register_resources(self, calculation_id: str, parameters: dict) -> None:
        """Register calculation resources with the resource manager."""
        if self.resource_manager is None:
            return

        user_cpu_cores = parameters.get('cpu_cores') or 1
        user_memory_mb = parameters.get('memory_mb') or get_memory_for_method(
            parameters.get('calculation_method', 'DFT'))
        calculation_method = parameters.get('calculation_method', 'DFT')

        try:
            self.resource_manager.register_calculation(
                calculation_id=calculation_id,
                cpu_cores=user_cpu_cores,
                memory_mb=user_memory_mb,
                calculation_method=calculation_method
            )
        except Exception as e:
            logger.warning(f"Failed to register calculation resources: {e}")

    def unregister_resources(self, calculation_id: str) -> None:
        """Unregister calculation resources from the resource manager."""
        if self.resource_manager is None:
            return
        try:
            self.resource_manager.unregister_calculation(calculation_id)
        except Exception as e:
            logger.warning(f"Failed to unregister calculation resources: {e}")

    def process_queue(self, active_futures: Dict[str, Future],
                      start_callback: Callable, error_callback: Callable,
                      is_shutdown: bool) -> None:
        """
        Process the calculation queue, starting calculations when resources allow.

        Args:
            active_futures: Dict of currently active calculation futures
            start_callback: Called with (QueuedCalculation) to start a calculation.
                           Should return True on success, False on failure.
            error_callback: Called with (calculation_id, error_message) on errors.
            is_shutdown: Whether the process manager is shutting down
        """
        with self._queue_lock:
            if self._queue_processing:
                logger.debug("Queue processing already in progress, skipping")
                return
            self._queue_processing = True

        try:
            logger.info(f"Processing calculation queue - Active: {len(active_futures)}/{self.max_parallel_instances}, Queued: {len(self.calculation_queue)}")
            started_count = 0

            while self.calculation_queue and not is_shutdown:
                started = False

                for i, queued_calc in enumerate(self.calculation_queue):
                    allocation_status, reason = self.check_and_decide(
                        queued_calc.parameters,
                        len(active_futures),
                        self.max_parallel_instances
                    )

                    if allocation_status == AllocationStatus.INSUFFICIENT_RESOURCES:
                        next_calc = self.calculation_queue.pop(i)
                        logger.warning(f"Removing calculation {next_calc.calculation_id} from queue due to insufficient resources: {reason}")
                        error_callback(next_calc.calculation_id, reason)
                        started = True
                        break

                    elif allocation_status == AllocationStatus.CAN_START:
                        next_calc = self.calculation_queue.pop(i)

                        self.register_resources(next_calc.calculation_id, next_calc.parameters)

                        if is_shutdown:
                            logger.warning(f"Executor has been shut down. Cannot schedule calculation {next_calc.calculation_id}")
                            self.unregister_resources(next_calc.calculation_id)
                            error_callback(next_calc.calculation_id, 'Process manager has been shut down')
                            started = True
                            break

                        success = start_callback(next_calc)
                        if success:
                            started_count += 1
                            user_cpu_cores = next_calc.parameters.get('cpu_cores') or 1
                            user_memory_mb = next_calc.parameters.get('memory_mb') or get_memory_for_method(
                                next_calc.parameters.get('calculation_method', 'DFT'))
                            logger.info(f"Started queued calculation {next_calc.calculation_id} with {user_cpu_cores} CPU cores and {user_memory_mb} MB memory ({len(active_futures)}/{self.max_parallel_instances} slots used, {len(self.calculation_queue)} remaining in queue)")
                        else:
                            self.unregister_resources(next_calc.calculation_id)
                            error_callback(next_calc.calculation_id, 'Failed to start calculation from queue')
                        started = True
                        break

                    else:
                        logger.debug(f"Calculation {queued_calc.calculation_id} should remain in queue: {reason}")

                if not started:
                    if self.calculation_queue:
                        logger.info(f"Queue processing paused - {len(self.calculation_queue)} calculations waiting for resources")
                        for j, queued_calc in enumerate(self.calculation_queue[:3]):
                            logger.debug(f"  Waiting calculation {j+1}: {queued_calc.calculation_id} (reason: {queued_calc.waiting_reason})")
                    break

        finally:
            with self._queue_lock:
                self._queue_processing = False

            if started_count > 0:
                logger.info(f"Queue processing completed - started {started_count} calculations from queue")

    def set_max_parallel_instances(self, max_instances: int, max_workers: int) -> int:
        """Update max parallel instances. Returns the new value."""
        old_max = self.max_parallel_instances
        self.max_parallel_instances = max(1, min(max_instances, max_workers))
        logger.info(f"Updated max parallel instances from {old_max} to {self.max_parallel_instances}")
        return self.max_parallel_instances

    def get_queued_calculations(self) -> list:
        """Get list of queued calculation IDs."""
        return [calc.calculation_id for calc in self.calculation_queue]

    def get_queue_status(self, active_count: int) -> dict:
        """Get current queue status information."""
        return {
            'active_calculations': active_count,
            'queued_calculations': len(self.calculation_queue),
            'max_parallel_instances': self.max_parallel_instances,
        }

    @property
    def has_queued(self) -> bool:
        """Check if there are calculations in the queue."""
        return bool(self.calculation_queue)
