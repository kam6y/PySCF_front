class TestSettingsAPI:
    def _settings_payload(self, calculations_directory='/tmp/pyscf-test'):
        return {
            'max_parallel_instances': 4,
            'max_cpu_utilization_percent': 95.0,
            'max_memory_utilization_percent': 95.0,
            'gpu_acceleration_enabled': False,
            'system_total_cores': 8,
            'system_total_memory_mb': 16384,
            'calculations_directory': calculations_directory,
            'timezone': 'UTC',
            'gemini_api_key': None,
            'research_email': None,
        }

    def test_get_settings_success(self, client, mocker):
        mock_settings = self._settings_payload()
        mock_service = mocker.patch('api.settings.get_settings_service')
        mock_service.return_value.get_settings.return_value = mock_settings

        response = client.get('/api/settings')

        assert response.status_code == 200
        assert response.json()['success'] is True
        assert response.json()['data']['settings'] == mock_settings

    def test_update_settings_success(self, client, mocker):
        updated = self._settings_payload('/tmp/pyscf-updated')
        mock_service = mocker.patch('api.settings.get_settings_service')
        mock_service.return_value.update_settings.return_value = updated

        response = client.put(
            '/api/settings',
            json=self._settings_payload('/tmp/pyscf-updated'),
        )

        assert response.status_code == 200
        assert response.json()['success'] is True
        assert response.json()['data']['settings'] == updated
        mock_service.return_value.update_settings.assert_called_once()
