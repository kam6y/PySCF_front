"""
Unit tests for SystemService diagnostics.
"""

import services.system_service as system_service
from generated_models import AppSettings
from services.system_service import SystemService


def test_get_system_diagnostics_masks_sensitive_settings(tmp_path, monkeypatch):
    """
    GIVEN current settings contain sensitive values
    WHEN system diagnostics are generated
    THEN settings diagnostics do not expose the raw sensitive values
    """
    sensitive_api_key = "plain-gemini-api-key"
    sensitive_email = "researcher@example.com"
    settings = AppSettings(
        max_parallel_instances=1,
        max_cpu_utilization_percent=95.0,
        max_memory_utilization_percent=95.0,
        system_total_cores=1,
        system_total_memory_mb=1024,
        calculations_directory=str(tmp_path),
        timezone="UTC",
        gemini_api_key=sensitive_api_key,
        research_email=sensitive_email,
        gpu_acceleration_enabled=False,
    )

    monkeypatch.setattr(system_service, "get_current_settings", lambda: settings)

    diagnostics = SystemService().get_system_diagnostics()
    diagnostics_settings = diagnostics["settings"]["settings"]

    assert diagnostics_settings["gemini_api_key"] == "***"
    assert diagnostics_settings["research_email"] == "***"
    assert sensitive_api_key not in str(diagnostics_settings)
    assert sensitive_email not in str(diagnostics_settings)
