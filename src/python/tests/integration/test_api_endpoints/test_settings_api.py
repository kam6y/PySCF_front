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
