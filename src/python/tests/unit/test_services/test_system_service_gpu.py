"""
Unit tests for SystemService GPU-related helpers.
"""

import importlib
import site
import subprocess
import sys
from unittest.mock import MagicMock

import pytest

import services.system_service as system_service
from services.exceptions import ValidationError
from services.system_service import SystemService


def test_is_module_available_adds_user_site_path(monkeypatch):
    """
    GIVEN user site-packages enabled and not in sys.path
    WHEN _is_module_available is called
    THEN user site-packages is appended to sys.path
    """
    service = SystemService()
    user_site = "/tmp/user-site"

    original_path = list(sys.path)
    if user_site in sys.path:
        sys.path.remove(user_site)

    try:
        monkeypatch.setattr(site, "ENABLE_USER_SITE", True)
        monkeypatch.setattr(site, "getusersitepackages", lambda: user_site)
        monkeypatch.setattr(importlib, "invalidate_caches", lambda: None)
        monkeypatch.setattr(system_service.util, "find_spec", lambda _: None)

        assert service._is_module_available("gpu4pyscf") is False
        assert user_site in sys.path
    finally:
        sys.path[:] = original_path


def test_get_gpu4pyscf_status_nvcc_missing(monkeypatch):
    """
    GIVEN Linux and nvcc is missing
    WHEN get_gpu4pyscf_status is called
    THEN CUDA detection is false and unsupported is reported
    """
    service = SystemService()

    monkeypatch.setattr(system_service.sys, "platform", "linux")
    monkeypatch.setattr(
        service,
        "_detect_cuda_version",
        lambda: (None, None, None, "nvcc command not found"),
    )
    monkeypatch.setattr(service, "_is_module_available", lambda _: False)

    status = service.get_gpu4pyscf_status()

    assert status["cuda_detected"] is False
    assert status["cuda_supported"] is False
    assert status["cuda_detection_message"] == "nvcc command not found"


def test_install_gpu4pyscf_with_confirm_false_raises_validation_error():
    """
    GIVEN confirm_install is False (default)
    WHEN install_gpu4pyscf is called
    THEN ValidationError is raised with the confirmation prompt message
    """
    # Arrange
    service = SystemService()

    # Act & Assert
    with pytest.raises(ValidationError) as exc_info:
        service.install_gpu4pyscf(confirm_install=False)

    assert "confirm_install=true" in str(exc_info.value)
    assert "mutates the Python environment" in str(exc_info.value)


def test_install_gpu4pyscf_with_confirm_true_proceeds_to_install(monkeypatch):
    """
    GIVEN confirm_install is True and Linux with supported CUDA
    WHEN install_gpu4pyscf is called
    THEN it proceeds past the confirmation gate and attempts installation
         (subprocess.run is mocked so no real pip install occurs)
    """
    # Arrange
    service = SystemService()

    monkeypatch.setattr(system_service.sys, "platform", "linux")
    monkeypatch.setattr(
        service,
        "_detect_cuda_version",
        lambda: ("12.4", 12, 4, None),
    )
    monkeypatch.setattr(service, "_is_site_writable", lambda: True)

    mock_run = MagicMock(
        return_value=subprocess.CompletedProcess(
            args=[], returncode=0, stdout="Successfully installed", stderr=""
        )
    )
    monkeypatch.setattr(system_service.subprocess, "run", mock_run)

    monkeypatch.setattr(service, "_is_module_available", lambda _: True)
    monkeypatch.setattr(service, "_get_distribution_version", lambda _names: "0.6.1")

    # Act
    result = service.install_gpu4pyscf(confirm_install=True)

    # Assert — reached the install path (subprocess.run was called)
    assert mock_run.call_count >= 1
    assert result["status"]["gpu4pyscf_installed"] is True
    # Verify pip command was built correctly (no shell=True)
    first_call_args = mock_run.call_args_list[0]
    assert isinstance(
        first_call_args[0][0], list
    ), "pip command must be a list, not a shell string"
