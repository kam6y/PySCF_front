"""Unit tests for resource manager fallback behavior."""

from unittest.mock import Mock

import quantum_calc.resource_manager as resource_manager
from quantum_calc.resource_manager import SystemResourceManager


def test_resource_manager_handles_unknown_logical_cpu_count(monkeypatch):
    """
    GIVEN psutil cannot determine logical CPU count
    WHEN SystemResourceManager initializes
    THEN resource constraints still contain a usable positive core count
    """
    mock_memory = Mock()
    mock_memory.total = 16 * 1024 * 1024 * 1024
    mock_memory.available = 8 * 1024 * 1024 * 1024
    mock_memory.percent = 50.0

    monkeypatch.setattr(resource_manager.psutil, "cpu_count", lambda logical=True: None)
    monkeypatch.setattr(resource_manager.psutil, "virtual_memory", lambda: mock_memory)
    monkeypatch.setattr(resource_manager.multiprocessing, "cpu_count", lambda: 3)

    manager = SystemResourceManager()

    constraints = manager.get_resource_constraints()
    assert constraints.system_total_cores == 3

    can_allocate, reason = manager.can_allocate_resources(
        cpu_cores=1,
        memory_mb=512,
        calculation_method="HF"
    )
    assert can_allocate is True
    assert "Resources available" in reason
