"""Unit tests for worker runtime calculator setup."""

from __future__ import annotations

from typing import Any
from unittest.mock import Mock

import pytest

from quantum_calc._worker_runtime import (
    _create_calculator_instance,
    _prepare_setup_parameters,
)


SUPPORTED_METHODS = ["DFT", "HF", "MP2", "CCSD", "CCSD_T", "TDDFT"]


def _build_recording_calculator() -> tuple[type, list[Any]]:
    instances: list[Any] = []

    class RecordingCalculator:
        def __init__(self, **kwargs: Any) -> None:
            self.kwargs = kwargs
            instances.append(self)

    return RecordingCalculator, instances


def _calculator_classes(calculator_class: type) -> dict[str, type]:
    return {method: calculator_class for method in SUPPORTED_METHODS}


@pytest.mark.parametrize(
    ("method", "expected_optimize_geometry"),
    [
        ("DFT", True),
        ("HF", True),
        ("MP2", True),
        ("CCSD", False),
        ("CCSD_T", False),
        ("TDDFT", False),
    ],
)
def test_create_calculator_instance_uses_method_default_optimize_geometry(
    method: str,
    expected_optimize_geometry: bool,
) -> None:
    """
    GIVEN no optimize_geometry parameter was supplied
    WHEN the worker creates a calculator instance
    THEN geometry optimization defaults only for methods that support it
    """
    calculator_class, instances = _build_recording_calculator()

    _create_calculator_instance(
        method,
        {"name": "Water", "calculation_method": method},
        "/tmp/calc",
        _calculator_classes(calculator_class),
        Mock(),
    )

    assert instances[0].kwargs["optimize_geometry"] is expected_optimize_geometry


def test_create_calculator_instance_preserves_explicit_optimize_geometry_false() -> None:
    """
    GIVEN a geometry-optimization method explicitly disables optimization
    WHEN the worker creates a calculator instance
    THEN the explicit value is preserved over the method default
    """
    calculator_class, instances = _build_recording_calculator()

    _create_calculator_instance(
        "DFT",
        {
            "name": "Water",
            "calculation_method": "DFT",
            "optimize_geometry": False,
        },
        "/tmp/calc",
        _calculator_classes(calculator_class),
        Mock(),
    )

    assert instances[0].kwargs["optimize_geometry"] is False


def test_create_calculator_instance_rejects_unknown_method() -> None:
    """
    GIVEN an unsupported calculation method reaches the worker
    WHEN the worker creates a calculator instance
    THEN it should fail instead of silently running DFT
    """
    calculator_class, _ = _build_recording_calculator()

    with pytest.raises(ValueError, match="Unsupported calculation method"):
        _create_calculator_instance(
            "UNKNOWN",
            {"name": "Water", "calculation_method": "UNKNOWN"},
            "/tmp/calc",
            _calculator_classes(calculator_class),
            Mock(),
        )


def test_prepare_setup_parameters_preserves_gpu_setting_snapshot() -> None:
    """
    GIVEN persisted job parameters include a GPU acceleration snapshot
    WHEN worker setup parameters are prepared
    THEN the calculator receives the persisted GPU setting
    """
    setup_params = _prepare_setup_parameters(
        {
            "calculation_method": "HF",
            "basis_function": "sto-3g",
            "charges": 0,
            "spin": 0,
            "solvent_method": "none",
            "solvent": "-",
            "gpu_acceleration_enabled": True,
        },
        memory_mb=1024,
    )

    assert setup_params["gpu_acceleration_enabled"] is True
