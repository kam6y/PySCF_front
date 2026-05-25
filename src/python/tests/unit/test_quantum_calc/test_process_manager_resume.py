"""
Unit tests for CalculationProcessManager pause/resume persistence.
"""

from types import SimpleNamespace

import pytest

from quantum_calc._status_transition import CalculationStatus
from quantum_calc._calculation_repository import CalculationRepository
from quantum_calc.process_manager import CalculationProcessManager


class FakeScheduler:
    def __init__(self):
        self.unregistered_calculation_ids = []
        self.processed = False

    def unregister_resources(self, calculation_id):
        self.unregistered_calculation_ids.append(calculation_id)

    def process_queue(self, **_kwargs):
        self.processed = True


class FakeStatusManager:
    def __init__(self):
        self.transitions = []
        self.notifications = []

    def transition(self, calculation_id, status, error_message=None):
        self.transitions.append((calculation_id, status, error_message))

    def notify(self, calculation_id, status, error_message=None):
        self.notifications.append((calculation_id, status, error_message))


class FailedFuture:
    def __init__(self, exception):
        self._exception = exception

    def exception(self):
        return self._exception


class RunningFuture:
    def __init__(self):
        self.cancel_called = False

    def done(self):
        return False

    def cancel(self):
        self.cancel_called = True
        return False


class FakeWorkerProcess:
    def __init__(self):
        self.terminate_called = False
        self.kill_called = False
        self.join_calls = []
        self._alive = True

    def is_alive(self):
        return self._alive

    def terminate(self):
        self.terminate_called = True
        self._alive = False

    def kill(self):
        self.kill_called = True
        self._alive = False

    def join(self, timeout=None):
        self.join_calls.append(timeout)


class FakeProcessPoolExecutor:
    def __init__(self, process):
        self._processes = {123: process}
        self.shutdown_calls = []

    def shutdown(self, wait=True, *, cancel_futures=False):
        self.shutdown_calls.append((wait, cancel_futures))


def test_pause_calculation_requires_active_future(tmp_path, monkeypatch):
    """
    GIVEN a persisted running calculation is not owned by an active future
    WHEN the process manager is asked to pause it directly
    THEN no pause flag or pausing transition is written
    """
    repository = CalculationRepository(base_dir=str(tmp_path))
    calc_id = "stale-running-calc"
    calc_dir = tmp_path / calc_id
    calc_dir.mkdir()
    repository.save_calculation_status(str(calc_dir), "running")

    manager = object.__new__(CalculationProcessManager)
    manager.active_futures = {}
    manager.status_manager = FakeStatusManager()
    manager._shutdown = True

    monkeypatch.setattr(
        "quantum_calc.get_current_settings",
        lambda: SimpleNamespace(calculations_directory=str(tmp_path)),
    )

    with pytest.raises(ValueError, match="no active worker"):
        manager.pause_calculation(calc_id)

    assert repository.read_calculation_status_details(str(calc_dir)) == (
        "running",
        None,
    )
    assert not (calc_dir / ".pause_requested").exists()
    assert manager.status_manager.transitions == []


def test_resume_calculation_persists_waiting_status(tmp_path, monkeypatch):
    """
    GIVEN a paused calculation is resumed but queued
    WHEN resume_calculation receives a waiting submit result
    THEN status.json is updated to waiting with the queue reason
    """
    repository = CalculationRepository(base_dir=str(tmp_path))
    calc_id = "paused-calc"
    calc_dir = tmp_path / calc_id
    calc_dir.mkdir()
    repository.save_calculation_parameters(
        str(calc_dir),
        {
            "name": "Paused Calc",
            "created_at": "2026-05-20T00:00:00",
            "calculation_method": "HF",
        },
    )
    repository.save_calculation_status(str(calc_dir), "paused")

    manager = object.__new__(CalculationProcessManager)
    manager._shutdown = True
    monkeypatch.setattr(
        "quantum_calc.get_current_settings",
        lambda: SimpleNamespace(calculations_directory=str(tmp_path)),
    )
    monkeypatch.setattr(
        manager,
        "submit_calculation",
        lambda calculation_id, params: (True, "waiting", "All slots are busy"),
    )

    result = manager.resume_calculation(calc_id)

    assert result == {
        "calculation_id": calc_id,
        "status": "waiting",
        "waiting_reason": "All slots are busy",
    }
    assert repository.read_calculation_status_details(str(calc_dir)) == (
        "waiting",
        "All slots are busy",
    )


def test_cleanup_future_with_worker_exception_persists_error_transition():
    """
    GIVEN a worker future failed before the worker wrote its result
    WHEN the parent process cleans up the future
    THEN the calculation status is transitioned to error instead of notification-only
    """
    calc_id = "crashed-calc"
    scheduler = FakeScheduler()
    status_manager = FakeStatusManager()

    manager = object.__new__(CalculationProcessManager)
    manager.active_futures = {calc_id: object()}
    manager.scheduler = scheduler
    manager.status_manager = status_manager
    manager.executor = None
    manager._resource_monitor_thread = None
    manager._shutdown = False

    manager._cleanup_future(calc_id, FailedFuture(RuntimeError("boom")))

    assert status_manager.transitions == [
        (calc_id, CalculationStatus.ERROR, "boom"),
    ]
    assert status_manager.notifications == []
    assert scheduler.unregistered_calculation_ids == [calc_id]
    assert scheduler.processed is True


def test_shutdown_force_terminates_running_worker_without_waiting():
    """
    GIVEN a running calculation future during application shutdown
    WHEN process manager shutdown is forced without waiting
    THEN worker processes are terminated and executor shutdown does not block.
    """
    future = RunningFuture()
    process = FakeWorkerProcess()
    executor = FakeProcessPoolExecutor(process)

    manager = object.__new__(CalculationProcessManager)
    manager.active_futures = {"running-calc": future}
    manager.executor = executor
    manager._resource_monitor_thread = None
    manager._shutdown = False

    manager.shutdown(wait=False, timeout=0.1, force=True)

    assert future.cancel_called is True
    assert process.terminate_called is True
    assert process.kill_called is False
    assert executor.shutdown_calls == [(False, True)]
    assert manager.executor is None
    assert manager.active_futures == {}
    assert manager._shutdown is True


def test_shutdown_process_manager_forwards_force_options(monkeypatch):
    """
    GIVEN the global process manager is initialized
    WHEN shutdown_process_manager receives forced shutdown options
    THEN they are forwarded to the concrete manager.
    """
    import quantum_calc.process_manager as process_manager_module

    class FakeManager:
        def __init__(self):
            self.shutdown_calls = []

        def shutdown(self, wait=True, timeout=None, force=False):
            self.shutdown_calls.append((wait, timeout, force))

    fake_manager = FakeManager()
    monkeypatch.setattr(process_manager_module, "_process_manager", fake_manager)

    process_manager_module.shutdown_process_manager(
        wait=False,
        timeout=0.2,
        force=True,
    )

    assert fake_manager.shutdown_calls == [(False, 0.2, True)]
    assert process_manager_module._process_manager is None
