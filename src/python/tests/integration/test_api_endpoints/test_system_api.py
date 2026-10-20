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
            ("/api/debug/system-diagnostics", "not enabled"),
            ("/api/debug/process-manager-diagnostics", "not enabled"),
            ("/api/debug/resource-manager-diagnostics", "not enabled"),
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


    @pytest.mark.parametrize(
        ("env_value", "env_action"),
        [
            (None, "del"),
            ("", "set"),
            ("prod", "set"),
        ],
        ids=[
            "env-absent",
            "env-empty",
            "env-typo-prod",
        ],
    )
    def test_debug_endpoint_returns_403_when_env_unrecognised(
        self,
        client,
        monkeypatch,
        env_value,
        env_action,
    ):
        """
        GIVEN PYSCF_ENV is absent, empty, or a typo (fail-closed default)
        WHEN GET /api/debug/system-diagnostics is called
        THEN 403 is returned with 'not enabled' message
        """
        # Arrange — override conftest defaults to exercise fail-closed path
        if env_action == "del":
            monkeypatch.delenv("PYSCF_ENV", raising=False)
        else:
            monkeypatch.setenv("PYSCF_ENV", env_value)
        monkeypatch.delenv("PYSCF_ENABLE_DEBUG_ENDPOINTS", raising=False)

        # Act
        response = client.get("/api/debug/system-diagnostics")

        # Assert
        assert response.status_code == 403
        data = response.json()
        assert data["success"] is False
        assert "not enabled" in data["error"].lower()

    @pytest.mark.parametrize(
        ("flag_value", "flag_action"),
        [
            (None, "del"),
            ("false", "set"),
        ],
        ids=[
            "flag-absent",
            "flag-false",
        ],
    )
    def test_debug_endpoint_returns_403_when_flag_not_enabled(
        self,
        client,
        monkeypatch,
        flag_value,
        flag_action,
    ):
        """
        GIVEN PYSCF_ENV is 'development' but PYSCF_ENABLE_DEBUG_ENDPOINTS
              is absent or falsy (Gate-2 failure, fail-closed)
        WHEN GET /api/debug/system-diagnostics is called
        THEN 403 is returned with 'not enabled' message
        """
        # Arrange — env passes Gate 1 but flag fails Gate 2
        monkeypatch.setenv("PYSCF_ENV", "development")
        if flag_action == "del":
            monkeypatch.delenv("PYSCF_ENABLE_DEBUG_ENDPOINTS", raising=False)
        else:
            monkeypatch.setenv("PYSCF_ENABLE_DEBUG_ENDPOINTS", flag_value)

        # Act
        response = client.get("/api/debug/system-diagnostics")

        # Assert
        assert response.status_code == 403
        data = response.json()
        assert data["success"] is False
        assert "not enabled" in data["error"].lower()

class TestDebugEndpointDoubleGate:
    """Unit tests for the fail-closed double-gate in _debug_endpoints_allowed.

    Gate 1: PYSCF_ENV must be in {"development", "test"}.
    Gate 2: PYSCF_ENABLE_DEBUG_ENDPOINTS must be truthy ("true" or "1").
    Both gates must pass; any other combination returns False (fail-closed).
    """

    @pytest.mark.parametrize(
        ("env_value", "flag_value", "expected"),
        [
            # Gate 1 fails: env absent or not in allowed set
            (None, None, False),
            (None, "true", False),
            ("", None, False),
            ("", "true", False),
            ("production", None, False),
            ("production", "true", False),
            ("prod", None, False),
            ("prod", "true", False),
            # Gate 1 passes, gate 2 fails: flag absent or falsy
            ("development", None, False),
            ("development", "", False),
            ("development", "false", False),
            ("development", "0", False),
            ("test", None, False),
            ("test", "0", False),
            # Both gates pass
            ("development", "true", True),
            ("development", "1", True),
            ("development", "TRUE", True),
            ("development", "True", True),
            ("test", "true", True),
            ("test", "1", True),
            # Case-insensitive env
            ("Development", "true", True),
            ("TEST", "1", True),
        ],
        ids=[
            "env-absent-flag-absent",
            "env-absent-flag-true",
            "env-empty-flag-absent",
            "env-empty-flag-true",
            "env-production-flag-absent",
            "env-production-flag-true",
            "env-prod-typo-flag-absent",
            "env-prod-typo-flag-true",
            "env-development-flag-absent",
            "env-development-flag-empty",
            "env-development-flag-false",
            "env-development-flag-zero",
            "env-test-flag-absent",
            "env-test-flag-zero",
            "env-development-flag-true",
            "env-development-flag-one",
            "env-development-flag-TRUE",
            "env-development-flag-True",
            "env-test-flag-true",
            "env-test-flag-one",
            "env-Development-flag-true",
            "env-TEST-flag-one",
        ],
    )
    def test_debug_endpoints_allowed_double_gate(
        self,
        monkeypatch,
        env_value,
        flag_value,
        expected,
    ):
        """
        GIVEN specific PYSCF_ENV and PYSCF_ENABLE_DEBUG_ENDPOINTS values
        WHEN _debug_endpoints_allowed is called
        THEN it returns the expected boolean (fail-closed)
        """
        from api.system import _debug_endpoints_allowed

        if env_value is None:
            monkeypatch.delenv("PYSCF_ENV", raising=False)
        else:
            monkeypatch.setenv("PYSCF_ENV", env_value)

        if flag_value is None:
            monkeypatch.delenv("PYSCF_ENABLE_DEBUG_ENDPOINTS", raising=False)
        else:
            monkeypatch.setenv("PYSCF_ENABLE_DEBUG_ENDPOINTS", flag_value)

        assert _debug_endpoints_allowed() is expected

    @pytest.mark.parametrize(
        ("value", "expected"),
        [
            (None, False),
            ("", False),
            ("false", False),
            ("0", False),
            ("true", True),
            ("1", True),
            ("TRUE", True),
            ("True", True),
            (" true ", True),
            (" 1 ", True),
        ],
        ids=[
            "none",
            "empty",
            "false-str",
            "zero-str",
            "true-lower",
            "one-str",
            "true-upper",
            "true-title",
            "true-whitespace",
            "one-whitespace",
        ],
    )
    def test_parse_debug_flag(self, value, expected):
        """
        GIVEN a raw flag string (or None)
        WHEN _parse_debug_flag is called
        THEN it returns the correct truthy/falsy interpretation
        """
        from api.system import _parse_debug_flag

        assert _parse_debug_flag(value) is expected


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


class TestGpu4PyscfPackagedModeGate:
    """SEC-004: Runtime install is disabled in packaged (production) builds."""

    def test_install_denied_in_production_without_opt_in(self, monkeypatch):
        """
        GIVEN PYSCF_ENV=production and PYSCF_ALLOW_RUNTIME_INSTALL is unset
        WHEN install_gpu4pyscf is called
        THEN ValidationError is raised before any installation attempt
        """
        from services.system_service import SystemService
        from services.exceptions import ValidationError as SvcValidationError

        monkeypatch.setenv("PYSCF_ENV", "production")
        monkeypatch.delenv("PYSCF_ALLOW_RUNTIME_INSTALL", raising=False)

        service = SystemService()
        with pytest.raises(SvcValidationError, match="disabled in packaged"):
            service.install_gpu4pyscf(confirm_install=True)

    def test_install_allowed_in_production_with_opt_in(self, monkeypatch):
        """
        GIVEN PYSCF_ENV=production and PYSCF_ALLOW_RUNTIME_INSTALL=1
        WHEN _runtime_install_allowed is called
        THEN it returns True (gate passes)
        """
        from services.system_service import SystemService

        monkeypatch.setenv("PYSCF_ENV", "production")
        monkeypatch.setenv("PYSCF_ALLOW_RUNTIME_INSTALL", "1")

        assert SystemService._runtime_install_allowed() is True

    def test_install_succeeds_in_production_with_opt_in(self, monkeypatch, mocker):
        """
        GIVEN PYSCF_ENV=production and PYSCF_ALLOW_RUNTIME_INSTALL=1
        WHEN install_gpu4pyscf is called with confirm_install=True
        THEN execution passes the gate and reaches _install_dependency_first
        """
        from services.system_service import (
            CUDA_PACKAGE_MAP,
            SystemService,
        )

        monkeypatch.setenv("PYSCF_ENV", "production")
        monkeypatch.setenv("PYSCF_ALLOW_RUNTIME_INSTALL", "1")

        service = SystemService()

        # Simulate Linux + CUDA 12.4
        monkeypatch.setattr("sys.platform", "linux")
        mocker.patch.object(
            service,
            "_detect_cuda_version",
            return_value=("12.4", 12, 4, None),
        )
        mocker.patch.object(service, "_is_site_writable", return_value=True)

        mock_install = mocker.patch.object(
            service,
            "_install_dependency_first",
            return_value=(True, ["gpu4pyscf-cuda12x==1.7.1"], "", ""),
        )
        mocker.patch.object(
            service,
            "get_gpu4pyscf_status",
            return_value={"gpu4pyscf_installed": True},
        )

        service.install_gpu4pyscf(confirm_install=True)

        # Verify execution reached _install_dependency_first (past the gate)
        mock_install.assert_called_once()
        gpu4pyscf_spec = mock_install.call_args.args[0]
        expected_pkg = CUDA_PACKAGE_MAP[12][0]
        assert gpu4pyscf_spec.startswith(expected_pkg)

    def test_install_denied_in_production_without_opt_in_via_gate(self, monkeypatch):
        """
        GIVEN PYSCF_ENV=production and PYSCF_ALLOW_RUNTIME_INSTALL is unset
        WHEN _runtime_install_allowed is called
        THEN it returns False (gate blocks)
        """
        from services.system_service import SystemService

        monkeypatch.setenv("PYSCF_ENV", "production")
        monkeypatch.delenv("PYSCF_ALLOW_RUNTIME_INSTALL", raising=False)

        assert SystemService._runtime_install_allowed() is False

    def test_install_allowed_in_development(self, monkeypatch):
        """
        GIVEN PYSCF_ENV=development
        WHEN _runtime_install_allowed is called
        THEN it returns True (gate passes)
        """
        from services.system_service import SystemService

        monkeypatch.setenv("PYSCF_ENV", "development")
        monkeypatch.delenv("PYSCF_ALLOW_RUNTIME_INSTALL", raising=False)

        assert SystemService._runtime_install_allowed() is True

    def test_install_denied_when_env_unset(self, monkeypatch):
        """
        GIVEN PYSCF_ENV is unset (defaults to empty string)
        WHEN _runtime_install_allowed is called
        THEN it returns False (fail-closed: unknown environments are restrictive)
        """
        from services.system_service import SystemService

        monkeypatch.delenv("PYSCF_ENV", raising=False)
        monkeypatch.delenv("PYSCF_ALLOW_RUNTIME_INSTALL", raising=False)

        assert SystemService._runtime_install_allowed() is False

    def test_install_denied_when_env_typo(self, monkeypatch):
        """
        GIVEN PYSCF_ENV is set to a typo like 'prod'
        WHEN _runtime_install_allowed is called
        THEN it returns False (fail-closed: unknown values are restrictive)
        """
        from services.system_service import SystemService

        monkeypatch.setenv("PYSCF_ENV", "prod")
        monkeypatch.delenv("PYSCF_ALLOW_RUNTIME_INSTALL", raising=False)

        assert SystemService._runtime_install_allowed() is False


class TestGpu4PyscfVersionPinning:
    """SEC-004: Verify gpu4pyscf top-level package is version-pinned."""

    def test_gpu4pyscf_install_spec_contains_version_pin(self):
        """
        GIVEN CUDA_PACKAGE_MAP and GPU4PYSCF_VERSION_BY_CUDA are defined
        WHEN _build_gpu4pyscf_install_spec is called for each CUDA generation
        THEN every resulting install spec contains '==' (version-pinned).
        """
        from services.system_service import CUDA_PACKAGE_MAP, SystemService

        for cuda_major, (package_name, _) in CUDA_PACKAGE_MAP.items():
            spec = SystemService._build_gpu4pyscf_install_spec(package_name, cuda_major)
            assert (
                "==" in spec
            ), f"CUDA {cuda_major}: spec '{spec}' is not version-pinned"
            assert spec.startswith(package_name)

    def test_gpu4pyscf_install_spec_unpinned_fallback(self):
        """
        GIVEN a CUDA generation that has no entry in GPU4PYSCF_VERSION_BY_CUDA
        WHEN _build_gpu4pyscf_install_spec is called
        THEN the bare package name is returned without '=='.
        """
        from services.system_service import SystemService

        spec = SystemService._build_gpu4pyscf_install_spec("gpu4pyscf-cuda99x", 99)
        assert spec == "gpu4pyscf-cuda99x"
        assert "==" not in spec

    def test_install_gpu4pyscf_passes_pinned_spec_to_install_dependency_first(
        self, monkeypatch, mocker
    ):
        """
        GIVEN a Linux platform with CUDA 12 detected
        WHEN install_gpu4pyscf is called
        THEN _install_dependency_first receives the '==' pinned spec matching
             GPU4PYSCF_VERSION_BY_CUDA[12], not the bare package name.
        """
        from services.system_service import (
            CUDA_PACKAGE_MAP,
            GPU4PYSCF_VERSION_BY_CUDA,
            SystemService,
        )

        monkeypatch.setenv("PYSCF_ENV", "development")

        service = SystemService()

        # Simulate Linux + CUDA 12.4
        monkeypatch.setattr("sys.platform", "linux")
        mocker.patch.object(
            service,
            "_detect_cuda_version",
            return_value=("12.4", 12, 4, None),
        )
        mocker.patch.object(service, "_is_site_writable", return_value=True)

        # Mock _install_dependency_first to capture its arguments
        mock_install = mocker.patch.object(
            service,
            "_install_dependency_first",
            return_value=(True, ["gpu4pyscf-cuda12x==1.7.1"], "", ""),
        )
        mocker.patch.object(
            service,
            "get_gpu4pyscf_status",
            return_value={"gpu4pyscf_installed": True},
        )

        service.install_gpu4pyscf(confirm_install=True)

        # Assert the pinned spec was passed
        mock_install.assert_called_once()
        gpu4pyscf_spec_arg = mock_install.call_args.args[0]
        expected_pkg = CUDA_PACKAGE_MAP[12][0]
        expected_version = GPU4PYSCF_VERSION_BY_CUDA[12]
        assert (
            "==" in gpu4pyscf_spec_arg
        ), f"Spec '{gpu4pyscf_spec_arg}' is not version-pinned"
        assert (
            gpu4pyscf_spec_arg == f"{expected_pkg}=={expected_version}"
        ), f"Expected '{expected_pkg}=={expected_version}', got '{gpu4pyscf_spec_arg}'"
