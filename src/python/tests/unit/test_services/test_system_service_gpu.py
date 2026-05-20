"""
Unit tests for SystemService GPU-related helpers.
"""

import importlib
import site
import sys

import services.system_service as system_service
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
