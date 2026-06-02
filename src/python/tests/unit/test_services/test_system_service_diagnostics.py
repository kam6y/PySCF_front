"""
Unit tests for SystemService diagnostics.
"""


import pytest

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


class TestSanitizePipOutput:
    """Unit tests for SystemService._sanitize_pip_output static method."""

    # ---- (A) Path redaction: 6 cases, each asserts path removed + placeholder present

    @pytest.mark.parametrize(
        ("text", "raw_path"),
        [
            (
                "Installing to /Users/alice/project/venv/lib/python3.11",
                "/Users/alice/project/venv/lib/python3.11",
            ),
            (
                "Found at /home/bob/.local/lib/site-packages",
                "/home/bob/.local/lib/site-packages",
            ),
            (
                "Mounted volume /mnt/data/packages/gpu4pyscf",
                "/mnt/data/packages/gpu4pyscf",
            ),
            (
                "Workspace path /workspace/app/lib/python3.11",
                "/workspace/app/lib/python3.11",
            ),
            (
                "Installed in /.venv/lib/python3.11/site-packages",
                "/.venv/lib/python3.11/site-packages",
            ),
            (
                'File "/usr/lib/python3.11/subprocess.py", line 1, in <module>',
                "/usr/lib/python3.11/subprocess.py",
            ),
        ],
        ids=[
            "absolute-users",
            "absolute-home",
            "absolute-mnt",
            "absolute-workspace",
            "venv-path",
            "traceback-file-path",
        ],
    )
    def test_sanitize_pip_output_redacts_path(self, text, raw_path):
        """
        GIVEN pip output containing a filesystem path
        WHEN _sanitize_pip_output is called
        THEN the path is replaced with <redacted-path>
        """
        result = SystemService._sanitize_pip_output(text)

        assert (
            raw_path not in result
        ), f"Path '{raw_path}' should be redacted in: {result}"
        assert "<redacted-path>" in result

    # ---- Multiple-paths test stays separate (no <redacted-path> assertion)

    def test_sanitize_pip_output_multiple_paths_all_redacted(self):
        """
        GIVEN pip output containing multiple different paths
        WHEN _sanitize_pip_output is called
        THEN all paths are redacted
        """
        # Arrange
        text = (
            "Source: /home/user/src/gpu4pyscf "
            "Target: /opt/conda/lib/python3.11/site-packages"
        )

        # Act
        result = SystemService._sanitize_pip_output(text)

        # Assert
        assert "/home/user/src/gpu4pyscf" not in result
        assert "/opt/conda/lib/python3.11/site-packages" not in result

    # ---- URL preservation tests (kept as-is: different assertion sets)

    def test_sanitize_pip_output_preserves_https_urls(self):
        """
        GIVEN pip output containing HTTPS URLs alongside filesystem paths
        WHEN _sanitize_pip_output is called
        THEN the URLs are preserved intact
        """
        # Arrange
        url = "https://pypi.org/simple/gpu4pyscf/"
        text = f"Downloading from {url} to /home/user/cache/pip"

        # Act
        result = SystemService._sanitize_pip_output(text)

        # Assert
        assert url in result
        assert "/home/user/cache/pip" not in result

    def test_sanitize_pip_output_preserves_http_urls(self):
        """
        GIVEN pip output containing HTTP URLs
        WHEN _sanitize_pip_output is called
        THEN the HTTP URLs are preserved intact
        """
        # Arrange
        url = "http://mirror.example.com/pypi/packages/gpu4pyscf-0.6.1.tar.gz"
        text = f"Fetching {url}"

        # Act
        result = SystemService._sanitize_pip_output(text)

        # Assert
        assert url in result

    # ---- (B) Passthrough tests: empty + no-match

    @pytest.mark.parametrize(
        "text",
        [
            "",
            "Successfully installed gpu4pyscf-cuda12x-0.6.1",
        ],
        ids=["empty-input", "no-match"],
    )
    def test_sanitize_pip_output_passthrough(self, text):
        """
        GIVEN text with no filesystem paths to redact
        WHEN _sanitize_pip_output is called
        THEN the text is returned unchanged
        """
        result = SystemService._sanitize_pip_output(text)

        assert result == text
