"""
Integration tests for System API endpoints.

Covers GPU4PySCF status detection and installation endpoints.
"""


import pytest

from generated_models import AppSettings


class TestSystemDiagnosticsAPI:
    """Integration tests for /api/debug/system-diagnostics endpoint."""

    def test_get_system_diagnostics_masks_sensitive_settings(
        self,
        client,
        mocker,
        tmp_path,
    ):
        """
        GIVEN current settings contain sensitive values
        WHEN GET /api/debug/system-diagnostics is called
        THEN settings diagnostics mask those values in the response payload
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
        mocker.patch(
            "services.system_service.get_current_settings",
            return_value=settings,
        )

        response = client.get("/api/debug/system-diagnostics")

        assert response.status_code == 200
        response_text = response.text
        data = response.json()
        settings_payload = data["data"]["settings"]["settings"]
        assert data["success"] is True
        assert settings_payload["gemini_api_key"] == "***"
        assert settings_payload["research_email"] == "***"
        assert sensitive_api_key not in response_text
        assert sensitive_email not in response_text


class TestDebugEndpointProductionGating:
    """Integration tests for debug endpoint production-environment gating."""

    @pytest.mark.parametrize(
        ("path", "error_contains"),
        [
            ("/api/debug/system-diagnostics", "production"),
            ("/api/debug/process-manager-diagnostics", None),
            ("/api/debug/resource-manager-diagnostics", None),
        ],
        ids=[
            "system-diagnostics",
            "process-manager-diagnostics",
            "resource-manager-diagnostics",
        ],
    )
    def test_debug_endpoint_returns_403_in_production(
        self,
        client,
        monkeypatch,
        path,
        error_contains,
    ):
        """
        GIVEN PYSCF_ENV is set to 'production'
        WHEN GET debug endpoint is called with valid auth
        THEN 403 is returned with success=False

        NOTE: A valid PYSCF_AUTH_TOKEN + matching header is required so the
        auth middleware lets the request through to the endpoint, where the
        production-environment gate returns 403.
        """
        # Arrange
        auth_token = "test-production-token"
        monkeypatch.setenv("PYSCF_ENV", "production")
        monkeypatch.setenv("PYSCF_AUTH_TOKEN", auth_token)

        # Act
        response = client.get(path, headers={"X-Auth-Token": auth_token})

        # Assert
        assert response.status_code == 403
        data = response.json()
        assert data["success"] is False
        if error_contains is not None:
            assert error_contains in data["error"].lower()

    def test_debug_system_diagnostics_returns_200_in_development(
        self,
        client,
        mocker,
        tmp_path,
    ):
        """
        GIVEN PYSCF_ENV is set to 'development' (default in test fixture)
        WHEN GET /api/debug/system-diagnostics is called
        THEN 200 is returned with diagnostic data
        """
        # Arrange — the conftest app fixture already sets PYSCF_ENV=development
        settings = AppSettings(
            max_parallel_instances=1,
            max_cpu_utilization_percent=95.0,
            max_memory_utilization_percent=95.0,
            system_total_cores=1,
            system_total_memory_mb=1024,
            calculations_directory=str(tmp_path),
            timezone="UTC",
            gemini_api_key="test-key",
            research_email="test@example.com",
            gpu_acceleration_enabled=False,
        )
        mocker.patch(
            "services.system_service.get_current_settings",
            return_value=settings,
        )

        # Act
        response = client.get("/api/debug/system-diagnostics")

        # Assert
        assert response.status_code == 200
        data = response.json()
        assert data["success"] is True
        assert "data" in data


class TestGpu4PyscfStatusAPI:
    """Integration tests for /api/system/gpu4pyscf-status endpoint."""

    def test_get_gpu4pyscf_status_success(self, client, mocker):
        """
        GIVEN SystemService returns GPU4PySCF status
        WHEN GET /api/system/gpu4pyscf-status is called
        THEN 200 OK is returned with status payload
        """
        mock_status = {
            "is_linux": True,
            "cuda_detected": True,
            "cuda_version": "12.4",
            "cuda_major": 12,
            "cuda_minor": 4,
            "cuda_supported": True,
            "cuda_detection_message": None,
            "recommended_gpu4pyscf_package": "gpu4pyscf-cuda12x",
            "recommended_cutensor_package": "cutensor-cu12",
            "gpu4pyscf_installed": True,
            "gpu4pyscf_version": "0.6.1",
            "cutensor_installed": True,
            "cutensor_version": "2.2.0",
        }
        mock_service = mocker.patch("api.system.get_system_service")
        mock_service.return_value.get_gpu4pyscf_status.return_value = mock_status

        response = client.get("/api/system/gpu4pyscf-status")

        assert response.status_code == 200
        data = response.json()
        assert data["success"] is True
        assert data["data"]["cuda_supported"] is True
        assert data["data"]["gpu4pyscf_installed"] is True


class TestGpu4PyscfInstallAPI:
    """Integration tests for /api/system/gpu4pyscf-install endpoint."""

    def test_install_gpu4pyscf_local_success(self, client, mocker):
        """
        GIVEN a local request and SystemService installs GPU4PySCF
        WHEN POST /api/system/gpu4pyscf-install is called
        THEN 200 OK is returned with installation result
        """
        mock_result = {
            "status": {
                "is_linux": True,
                "cuda_detected": True,
                "cuda_supported": True,
                "gpu4pyscf_installed": True,
                "cutensor_installed": True,
            },
            "packages": ["gpu4pyscf-cuda12x", "cutensor-cu12"],
            "used_user_site": True,
            "pip_stdout": "",
            "pip_stderr": "",
        }
        mock_get_host = mocker.patch(
            "api.system._get_client_host",
            return_value="127.0.0.1",
        )
        mock_service = mocker.patch("api.system.get_system_service")
        mock_service.return_value.install_gpu4pyscf.return_value = mock_result

        response = client.post("/api/system/gpu4pyscf-install")

        mock_get_host.assert_called_once()
        assert response.status_code == 200
        data = response.json()
        assert data["success"] is True
        assert data["data"]["status"]["gpu4pyscf_installed"] is True
        mock_service.return_value.install_gpu4pyscf.assert_called_once_with(
            include_cutensor=True,
            force_reinstall=False,
            confirm_install=False,
        )

    def test_install_gpu4pyscf_blocks_remote(self, client, mocker):
        """
        GIVEN a non-local request
        WHEN POST /api/system/gpu4pyscf-install is called
        THEN 403 Forbidden is returned
        """
        mocker.patch("api.system._get_client_host", return_value="10.10.10.10")
        mock_service = mocker.patch("api.system.get_system_service")

        response = client.post("/api/system/gpu4pyscf-install")

        assert response.status_code == 403
        data = response.json()
        assert data["success"] is False
        mock_service.assert_not_called()
