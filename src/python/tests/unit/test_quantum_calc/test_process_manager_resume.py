"""
Unit tests for CalculationProcessManager pause/resume persistence.
"""

from types import SimpleNamespace

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
    callback_calls = []

    manager = object.__new__(CalculationProcessManager)
    manager.active_futures = {calc_id: object()}
    manager.completion_callbacks = {
        calc_id: [
            lambda calculation_id, success, error_message: callback_calls.append(
                (calculation_id, success, error_message)
            )
        ]
    }
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
    assert callback_calls == [(calc_id, False, "boom")]
    assert scheduler.unregistered_calculation_ids == [calc_id]
    assert scheduler.processed is True
