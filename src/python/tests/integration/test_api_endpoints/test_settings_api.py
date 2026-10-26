class TestSettingsAPI:
    def _settings_payload(self, calculations_directory="/tmp/pyscf-test"):
        return {
            "max_parallel_instances": 4,
            "max_cpu_utilization_percent": 95.0,
            "max_memory_utilization_percent": 95.0,
            "gpu_acceleration_enabled": False,
            "system_total_cores": 8,
            "system_total_memory_mb": 16384,
            "calculations_directory": calculations_directory,
            "timezone": "UTC",
            "gemini_api_key": None,
            "research_email": None,
        }

    def test_get_settings_success(self, client, mocker):
        mock_settings = self._settings_payload()
        mock_service = mocker.patch("api.settings.get_settings_service")
        mock_service.return_value.get_settings.return_value = mock_settings

        response = client.get("/api/settings")

        assert response.status_code == 200
        assert response.json()["success"] is True
        # GET response is masked: real key is stripped, configured flag is added
        expected_settings = {
            **mock_settings,
            "gemini_api_key": "",
            "gemini_api_key_configured": bool(mock_settings.get("gemini_api_key")),
        }
        assert response.json()["data"]["settings"] == expected_settings

    def test_update_settings_success(self, client, mocker):
        payload = self._settings_payload("/tmp/pyscf-updated")
        updated = self._settings_payload("/tmp/pyscf-updated")
        mock_service = mocker.patch("api.settings.get_settings_service")
        mock_service.return_value.update_settings.return_value = updated

        response = client.put(
            "/api/settings",
            json=payload,
        )

        assert response.status_code == 200
        assert response.json()["success"] is True
        # PUT response is masked: real key is stripped, configured flag is added
        expected_updated = {
            **updated,
            "gemini_api_key": "",
            "gemini_api_key_configured": bool(updated.get("gemini_api_key")),
        }
        assert response.json()["data"]["settings"] == expected_updated
        mock_service.return_value.update_settings.assert_called_once_with(payload)

    def test_get_settings_masks_configured_api_key(self, client, mocker):
        """
        GIVEN the service returns settings with a non-empty gemini_api_key
        WHEN GET /api/settings is called
        THEN the response has gemini_api_key == "" and gemini_api_key_configured == True
        """
        mock_settings = {
            **self._settings_payload(),
            "gemini_api_key": "real-secret-key",
        }
        mock_service = mocker.patch("api.settings.get_settings_service")
        mock_service.return_value.get_settings.return_value = mock_settings

        response = client.get("/api/settings")

        assert response.status_code == 200
        assert response.json()["success"] is True
        settings = response.json()["data"]["settings"]
        assert settings["gemini_api_key"] == ""
        assert settings["gemini_api_key_configured"] is True
        assert "real-secret-key" not in str(response.json())

    def test_put_settings_masks_configured_api_key_in_response(self, client, mocker):
        """
        GIVEN the service returns updated settings with a non-empty gemini_api_key
        WHEN PUT /api/settings is called
        THEN the response has gemini_api_key == "" and gemini_api_key_configured == True
        """
        payload = self._settings_payload("/tmp/pyscf-updated")
        updated = {
            **self._settings_payload("/tmp/pyscf-updated"),
            "gemini_api_key": "sk-new-secret",
        }
        mock_service = mocker.patch("api.settings.get_settings_service")
        mock_service.return_value.update_settings.return_value = updated

        response = client.put("/api/settings", json=payload)

        assert response.status_code == 200
        assert response.json()["success"] is True
        settings = response.json()["data"]["settings"]
        assert settings["gemini_api_key"] == ""
        assert settings["gemini_api_key_configured"] is True
        assert "sk-new-secret" not in str(response.json())

    def test_put_settings_succeeds_without_gemini_api_key(self, client, mocker):
        """
        GIVEN a valid settings payload that omits gemini_api_key entirely
              (as the frontend initial-setup handler does)
        WHEN PUT /api/settings is called
        THEN 200 is returned with success=True, because gemini_api_key is
             optional (default None) and extra='forbid' only rejects unknown
             fields, not missing optional ones.
        """
        payload = self._settings_payload()
        del payload["gemini_api_key"]

        # Service returns settings with gemini_api_key defaulted to None
        updated = {**payload, "gemini_api_key": None}
        mock_service = mocker.patch("api.settings.get_settings_service")
        mock_service.return_value.update_settings.return_value = updated

        response = client.put("/api/settings", json=payload)

        assert response.status_code == 200
        assert response.json()["success"] is True

    def test_put_settings_rejects_unknown_fields(self, client):
        """
        GIVEN a valid settings payload with an extra unknown field
        WHEN PUT /api/settings is called
        THEN 400 is returned because AppSettings has extra='forbid'
        """
        payload = {**self._settings_payload(), "unknown_field": "should be rejected"}

        response = client.put("/api/settings", json=payload)

        assert response.status_code == 400
        body = response.json()
        assert body["success"] is False
        assert "unknown_field" in body["error"]

    def test_put_settings_rejects_response_only_field_gemini_api_key_configured(
        self, client
    ):
        """
        GIVEN a settings payload containing 'gemini_api_key_configured' (a
              response-only field added by _mask_api_key_for_response, not a
              declared AppSettings field)
        WHEN PUT /api/settings is called
        THEN 400 is returned because AppSettings has extra='forbid'

        This documents the App.tsx coupling where the GET response includes
        gemini_api_key_configured but the request schema forbids it, and
        would catch a regression if the response-masking ever changes.
        """
        payload = {**self._settings_payload(), "gemini_api_key_configured": True}

        response = client.put("/api/settings", json=payload)

        assert response.status_code == 400
        body = response.json()
        assert body["success"] is False
        assert "gemini_api_key_configured" in body["error"]
