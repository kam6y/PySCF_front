"""Process pool manager for CPU-bound quantum chemistry calculations."""

import os
import logging
import multiprocessing
from concurrent.futures import ProcessPoolExecutor, Future
from typing import Dict, Any, Optional, Callable, List
import threading
import time

from .exceptions import ProcessManagerError, PauseRequestedException
from .pause_manager import pause_manager
from ._worker_runtime import calculation_worker, _worker_initializer
from ._queue_scheduler import CalculationQueueScheduler
from ._status_transition import CalculationStatusManager, CalculationStatus

logger = logging.getLogger(__name__)


class CalculationProcessManager:
    """Manages a process pool for quantum chemistry calculations with queueing support and resource management."""

    def __init__(self, max_workers: Optional[int] = None, max_parallel_instances: Optional[int] = None,
                 notification_callback: Optional[Callable] = None):
        self.max_workers = max_workers or multiprocessing.cpu_count()
        self.active_futures: Dict[str, Future] = {}
        self.completion_callbacks: Dict[str, List[Callable]] = {}
        self._shutdown = False

        # Initialize resource manager
        resource_manager = None
        try:
            from .resource_manager import get_resource_manager
            resource_manager = get_resource_manager()
            logger.info("Resource manager initialized successfully")
        except Exception as e:
            logger.error(f"Failed to initialize resource manager: {e}")
            logger.warning("Continuing without resource manager - calculations may not be resource-constrained")

        # Initialize status manager
        self.status_manager = CalculationStatusManager(notification_callback=notification_callback)

        # Initialize queue scheduler
        effective_max_parallel = max_parallel_instances or self.max_workers
        self.scheduler = CalculationQueueScheduler(
            resource_manager=resource_manager,
            max_parallel_instances=effective_max_parallel
        )

        # Background resource monitoring
        self._resource_monitor_thread = None
        self._resource_monitor_stop_event = threading.Event()
        self._resource_monitor_interval = 10.0

        # Initialize process pool executor
        self.executor = None
        try:
            self.executor = ProcessPoolExecutor(
                max_workers=self.max_workers,
                mp_context=multiprocessing.get_context('spawn'),
                initializer=_worker_initializer
            )
            logger.info(f"ProcessPoolExecutor created successfully with {self.max_workers} workers")
        except Exception as e:
            logger.error(f"Failed to create ProcessPoolExecutor: {e}")
            raise ProcessManagerError(f"Cannot create process pool: {str(e)}")

        if resource_manager is not None:
            try:
                self._start_resource_monitoring()
                logger.info("Background resource monitoring started")
            except Exception as e:
                logger.warning(f"Failed to start resource monitoring: {e}")
        else:
            logger.warning("Skipping resource monitoring due to resource manager initialization failure")

        logger.info(f"Initialized CalculationProcessManager with {self.max_workers} worker processes and {self.scheduler.max_parallel_instances} max parallel instances")

    # --- Properties for backward compatibility ---

    @property
    def notification_callback(self):
        return self.status_manager.notification_callback

    @notification_callback.setter
    def notification_callback(self, value):
        self.status_manager.notification_callback = value

    @property
    def max_parallel_instances(self):
        return self.scheduler.max_parallel_instances

    @max_parallel_instances.setter
    def max_parallel_instances(self, value):
        self.scheduler.max_parallel_instances = value

    @property
    def resource_manager(self):
        return self.scheduler.resource_manager

    @property
    def calculation_queue(self):
        return self.scheduler.calculation_queue

    # --- Resource Monitoring ---

    def _start_resource_monitoring(self):
        if self._resource_monitor_thread is None or not self._resource_monitor_thread.is_alive():
            self._resource_monitor_stop_event.clear()
            self._resource_monitor_thread = threading.Thread(
                target=self._resource_monitoring_loop,
                name="ResourceMonitor",
                daemon=True
            )
            self._resource_monitor_thread.start()
            logger.info("Started background resource monitoring")

    def _stop_resource_monitoring(self):
        if self._resource_monitor_thread and self._resource_monitor_thread.is_alive():
            logger.info("Stopping background resource monitoring")
            self._resource_monitor_stop_event.set()
            self._resource_monitor_thread.join(timeout=5.0)
            if self._resource_monitor_thread.is_alive():
                logger.warning("Resource monitoring thread did not stop cleanly")

    def _resource_monitoring_loop(self):
        logger.info("Resource monitoring loop started")
        while not self._resource_monitor_stop_event.is_set():
            try:
                if self.scheduler.has_queued and not self._shutdown:
                    should_process = False
                    if self.scheduler.resource_manager is not None:
                        try:
                            has_improved, reason = self.scheduler.resource_manager.has_resources_improved()
                            if has_improved:
                                logger.info(f"Resources improved: {reason}. Processing queue...")
                                should_process = True
                        except Exception as e:
                            logger.warning(f"Failed to check resource improvements: {e}. Processing queue anyway.")
                            should_process = True
                    else:
                        should_process = True
                    if should_process:
                        self._process_queue()
                self._resource_monitor_stop_event.wait(timeout=self._resource_monitor_interval)
            except Exception as e:
                logger.error(f"Error in resource monitoring loop: {e}")
                self._resource_monitor_stop_event.wait(timeout=self._resource_monitor_interval)
        logger.info("Resource monitoring loop stopped")

    # --- Core API ---

    def set_max_parallel_instances(self, max_instances: int) -> None:
        old_max = self.scheduler.max_parallel_instances
        self.scheduler.set_max_parallel_instances(max_instances, self.max_workers)
        if self.scheduler.max_parallel_instances > old_max:
            logger.info(f"Max parallel instances increased, processing queue (current queue size: {len(self.scheduler.calculation_queue)})")
            self._process_queue()

    def submit_calculation(self, calculation_id: str, parameters: dict) -> tuple[bool, str, Optional[str]]:
        if self._shutdown:
            logger.error("Cannot submit calculation: process manager is shut down")
            return False, 'error', None
        if self.executor is None:
            logger.error("Process pool executor is not available")
            return False, 'error', None

        allocation_status, reason = self.scheduler.check_and_decide(
            parameters, len(self.active_futures), self.scheduler.max_parallel_instances
        )

        if allocation_status == AllocationStatus.INSUFFICIENT_RESOURCES:
            logger.error(f"System resources insufficient for calculation {calculation_id}: {reason}")
            return False, 'error', reason

        if allocation_status == AllocationStatus.SHOULD_QUEUE:
            self.scheduler.enqueue(calculation_id, parameters, reason)
            self.status_manager.notify(calculation_id, 'waiting', None)
            return True, 'waiting', reason

        # CAN_START
        try:
            logger.info(f"About to submit calculation {calculation_id} to executor")
            self.scheduler.register_resources(calculation_id, parameters)

            future = self.executor.submit(calculation_worker, calculation_id, parameters)
            self.active_futures[calculation_id] = future
            future.add_done_callback(lambda f: self._cleanup_future(calculation_id, f))

            user_cpu_cores = parameters.get('cpu_cores') or 1
            user_memory_mb = parameters.get('memory_mb') or 0
            logger.info(f"Started calculation {calculation_id} immediately with {user_cpu_cores} CPU cores and {user_memory_mb} MB memory ({len(self.active_futures)}/{self.scheduler.max_parallel_instances} slots used)")

            self.status_manager.notify(calculation_id, 'running', None)
            return True, 'running', None

        except Exception as e:
            logger.error(f"Failed to submit calculation {calculation_id}: {e}")
            import traceback
            logger.error(f"Full traceback:\n{traceback.format_exc()}")
            self.scheduler.unregister_resources(calculation_id)
            return False, 'error', None

    def _cleanup_future(self, calculation_id: str, future: Future):
        success = False
        error_message = None

        try:
            if calculation_id in self.active_futures:
                del self.active_futures[calculation_id]

            self.scheduler.unregister_resources(calculation_id)

            if future.exception():
                exception = future.exception()
                if isinstance(exception, PauseRequestedException):
                    logger.info(f"Calculation {calculation_id} was paused by user request")
                    self.status_manager.transition(calculation_id, CalculationStatus.PAUSED)
                    pause_manager.clear_pause_request(calculation_id)
                    # Remove pause flag file
                    from quantum_calc._calculation_repository import CalculationRepository
                    from quantum_calc import get_current_settings
                    settings = get_current_settings()
                    repository = CalculationRepository(base_dir=settings.calculations_directory)
                    calc_dir = str(repository.resolve_calculation_path(calculation_id))
                    pause_manager.remove_pause_flag_file(calc_dir)
                else:
                    error_message = str(exception)
                    logger.error(f"Calculation {calculation_id} failed with exception: {error_message}")
                    self.status_manager.transition(
                        calculation_id,
                        CalculationStatus.ERROR,
                        error_message,
                    )
            else:
                success, calc_error = future.result()
                if success:
                    logger.info(f"Calculation {calculation_id} completed successfully")
                    self.status_manager.notify(calculation_id, 'completed', None)
                else:
                    error_message = calc_error
                    logger.warning(f"Calculation {calculation_id} failed: {error_message}")
                    self.status_manager.notify(calculation_id, 'error', error_message)

            # Call completion callbacks
            if calculation_id in self.completion_callbacks:
                callbacks = self.completion_callbacks.pop(calculation_id)
                for callback in callbacks:
                    try:
                        callback(calculation_id, success, error_message)
                    except Exception as callback_error:
                        logger.error(f"Error in completion callback for {calculation_id}: {callback_error}")

            logger.info(f"Calculation {calculation_id} cleanup completed. Processing queue for waiting calculations...")
            self._process_queue()

        except Exception as e:
            logger.error(f"Error cleaning up future for {calculation_id}: {e}")

    def _start_from_queue(self, queued_calc) -> bool:
        """Callback for queue scheduler to start a calculation from the queue.

        Returns True if executor.submit() succeeded (calculation is running),
        False if submission itself failed.
        """
        try:
            if self._shutdown:
                return False
            if self.executor is None:
                logger.error("Process pool executor is not available")
                return False

            future = self.executor.submit(
                calculation_worker, queued_calc.calculation_id, queued_calc.parameters
            )
            self.active_futures[queued_calc.calculation_id] = future
            future.add_done_callback(
                lambda f, calc_id=queued_calc.calculation_id: self._cleanup_future(calc_id, f)
            )
        except RuntimeError as e:
            if 'cannot schedule new futures after shutdown' in str(e):
                logger.error(f"Cannot schedule calculation {queued_calc.calculation_id}: executor has been shut down")
                return False
            raise
        except Exception as e:
            logger.error(f"Failed to start queued calculation {queued_calc.calculation_id}: {e}")
            return False

        # Calculation is now running — update status (best-effort; failure here
        # must not cause unregister since the worker is already executing)
        try:
            self.status_manager.transition(queued_calc.calculation_id, CalculationStatus.RUNNING)
        except Exception as e:
            logger.warning(f"Failed to update status for queued calculation {queued_calc.calculation_id}: {e}")

        return True

    def _handle_queue_error(self, calculation_id: str, error_message: str) -> None:
        """Callback for queue scheduler to handle errors."""
        self.status_manager.transition(calculation_id, CalculationStatus.ERROR, error_message)

    def _process_queue(self):
        self.scheduler.process_queue(
            active_futures=self.active_futures,
            start_callback=self._start_from_queue,
            error_callback=self._handle_queue_error,
            is_shutdown=self._shutdown
        )

    # --- Pause / Resume ---

    def pause_calculation(self, calculation_id: str) -> bool:
        logger.info(f"Pause requested for calculation: {calculation_id}")

        from quantum_calc._calculation_repository import CalculationRepository
        from quantum_calc import get_current_settings
        settings = get_current_settings()
        repository = CalculationRepository(base_dir=settings.calculations_directory)

        calc_dir = str(repository.resolve_calculation_path(calculation_id))
        if not os.path.exists(calc_dir):
            raise ValueError(f"Calculation not found: {calculation_id}")

        status, _ = repository.read_calculation_status_details(calc_dir)
        if status != 'running':
            raise ValueError(f"Calculation is not running (status: {status})")

        future = self.active_futures.get(calculation_id)
        if future is None or future.done():
            raise ValueError(
                f"Calculation has no active worker (calculation_id: {calculation_id})"
            )

        pause_manager.create_pause_flag_file(calc_dir)
        pause_manager.request_pause(calculation_id)

        self.status_manager.transition(calculation_id, CalculationStatus.PAUSING)

        logger.info(f"Pause request accepted for calculation: {calculation_id}")
        return True

    def resume_calculation(self, calculation_id: str) -> Dict[str, Any]:
        logger.info(f"Resuming calculation: {calculation_id}")

        from quantum_calc._calculation_repository import CalculationRepository
        from quantum_calc import get_current_settings
        settings = get_current_settings()
        repository = CalculationRepository(base_dir=settings.calculations_directory)

        calc_dir = str(repository.resolve_calculation_path(calculation_id))
        if not os.path.exists(calc_dir):
            raise ValueError(f"Calculation not found: {calculation_id}")

        status, _ = repository.read_calculation_status_details(calc_dir)
        if status != 'paused':
            raise ValueError(f"Calculation is not paused (status: {status})")

        pause_state = repository.load_pause_state(calc_dir)
        params = repository.read_calculation_parameters(calc_dir)
        if not params:
            raise ValueError("No calculation parameters found")

        params['resume_from_pause'] = True
        if pause_state:
            params['pause_state'] = pause_state

        success, message, waiting_reason = self.submit_calculation(calculation_id, params)
        if not success:
            raise ValueError(f"Failed to resume calculation: {message}")

        if message in {'waiting', 'running'}:
            repository.save_calculation_status(calc_dir, message, waiting_reason)

        logger.info(f"Calculation resumed: {calculation_id}")
        result = {'calculation_id': calculation_id, 'status': message}
        if waiting_reason is not None:
            result['waiting_reason'] = waiting_reason
        return result

    # --- Query methods ---

    def is_running(self, calculation_id: str) -> bool:
        future = self.active_futures.get(calculation_id)
        return future is not None and not future.done()

    def register_completion_callback(self, calculation_id: str, callback: Callable):
        if calculation_id not in self.completion_callbacks:
            self.completion_callbacks[calculation_id] = []
        self.completion_callbacks[calculation_id].append(callback)

    def get_active_calculations(self) -> list:
        return [calc_id for calc_id, future in self.active_futures.items() if not future.done()]

    def get_queued_calculations(self) -> list:
        return self.scheduler.get_queued_calculations()

    def get_queue_status(self) -> dict:
        status = self.scheduler.get_queue_status(len(self.active_futures))
        status['max_workers'] = self.max_workers
        return status

    def get_diagnostics(self) -> Dict[str, Any]:
        """Return a read-only diagnostic snapshot of internal state."""
        return {
            'is_shutdown': self._shutdown,
            'active_futures_count': len(self.active_futures),
            'active_calculation_ids': list(self.active_futures.keys()),
            'max_workers': self.max_workers,
            'max_parallel_instances': self.max_parallel_instances,
            'queued_calculations_count': len(self.calculation_queue),
            'completion_callbacks_count': len(self.completion_callbacks),
            'resource_monitoring': {
                'monitoring_active': (
                    self._resource_monitor_thread is not None
                    and self._resource_monitor_thread.is_alive()
                ),
                'monitoring_interval': self._resource_monitor_interval,
            },
            'executor_available': self.executor is not None,
            'executor_type': type(self.executor).__name__ if self.executor is not None else None,
        }

    # --- Lifecycle ---

    def _cancel_active_futures(self) -> None:
        for calculation_id, future in list(self.active_futures.items()):
            if future.done():
                continue
            if future.cancel():
                logger.info(f"Cancelled pending calculation future: {calculation_id}")
            else:
                logger.warning(
                    f"Calculation {calculation_id} is already running; terminating worker process"
                )

    def _terminate_executor_processes(self, timeout: float) -> None:
        if self.executor is None:
            return

        processes = getattr(self.executor, '_processes', None)
        if not processes:
            logger.info("No process pool workers to terminate")
            return

        worker_processes = [process for process in processes.values() if process is not None]
        if not worker_processes:
            logger.info("No process pool workers to terminate")
            return

        logger.warning(f"Terminating {len(worker_processes)} process pool workers")
        for process in worker_processes:
            if process.is_alive():
                process.terminate()

        deadline = time.monotonic() + timeout
        for process in worker_processes:
            remaining = max(0.0, deadline - time.monotonic())
            if process.is_alive():
                process.join(timeout=remaining)

        stubborn_processes = [process for process in worker_processes if process.is_alive()]
        if not stubborn_processes:
            return

        logger.warning(
            f"Killing {len(stubborn_processes)} process pool workers that ignored SIGTERM"
        )
        for process in stubborn_processes:
            kill = getattr(process, 'kill', None)
            if kill is not None:
                kill()
            else:
                process.terminate()
        for process in stubborn_processes:
            process.join(timeout=1.0)

    def shutdown(
        self,
        wait: bool = True,
        timeout: Optional[float] = None,
        force: bool = False,
    ) -> None:
        if self._shutdown:
            return
        self._shutdown = True
        self._stop_resource_monitoring()
        if self.executor:
            logger.info(
                f"Shutting down process pool with {len(self.active_futures)} active calculations"
            )
            try:
                if force:
                    self._cancel_active_futures()
                    self._terminate_executor_processes(timeout or 5.0)
                self.executor.shutdown(wait=wait, cancel_futures=force or not wait)
                logger.info("Process pool shut down successfully")
            except Exception as e:
                logger.error(f"Error during process pool shutdown: {e}")
            finally:
                self.executor = None
                self.active_futures.clear()

    def __del__(self):
        if not self._shutdown:
            self.shutdown(wait=False)


# --- Global process manager ---

from .resource_manager import AllocationStatus  # re-export for backward compat

_process_manager: Optional[CalculationProcessManager] = None


def initialize_process_manager_with_callback(notification_callback: Optional[Callable] = None,
                                             max_parallel_instances: Optional[int] = None) -> CalculationProcessManager:
    global _process_manager
    if _process_manager is not None:
        logger.warning("Process manager already initialized. Updating notification callback.")
        _process_manager.notification_callback = notification_callback
        if max_parallel_instances is not None:
            _process_manager.set_max_parallel_instances(max_parallel_instances)
        return _process_manager

    if max_parallel_instances is None:
        max_parallel_instances = min(4, multiprocessing.cpu_count())

    logger.info(f"Initializing process manager with notification callback, max_parallel={max_parallel_instances}")

    try:
        _process_manager = CalculationProcessManager(
            max_parallel_instances=max_parallel_instances,
            notification_callback=notification_callback
        )
        logger.info("Process manager initialized successfully with notification callback")
    except Exception as e:
        logger.error(f"Failed to initialize process manager: {e}")
        _process_manager = None
        raise ProcessManagerError(f"Failed to initialize process manager: {str(e)}")

    return _process_manager


def get_process_manager() -> CalculationProcessManager:
    global _process_manager
    if _process_manager is None:
        default_max_parallel = min(4, multiprocessing.cpu_count())
        logger.warning(f"Process manager not initialized with callback. Using default settings: max_parallel={default_max_parallel}")
        try:
            _process_manager = CalculationProcessManager(max_parallel_instances=default_max_parallel)
            logger.info("Process manager initialized with defaults (no callback)")
        except Exception as e:
            logger.error(f"Failed to initialize process manager: {e}")
            _process_manager = None
            raise ProcessManagerError(f"Failed to initialize process manager: {str(e)}")
    return _process_manager


def update_process_manager_settings():
    global _process_manager
    if _process_manager is not None:
        try:
            from .settings_manager import get_current_settings
            settings = get_current_settings()
            _process_manager.set_max_parallel_instances(settings.max_parallel_instances)
            logger.info(f"Updated process manager settings: max_parallel_instances={settings.max_parallel_instances}")
        except Exception as e:
            logger.warning(f"Failed to update process manager settings: {e}. Keeping current settings.")


def shutdown_process_manager(
    wait: bool = True,
    timeout: Optional[float] = None,
    force: bool = False,
) -> None:
    global _process_manager
    if _process_manager is not None:
        _process_manager.shutdown(wait=wait, timeout=timeout, force=force)
        _process_manager = None
