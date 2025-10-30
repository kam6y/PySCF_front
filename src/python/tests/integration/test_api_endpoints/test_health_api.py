"""
Integration tests for Health Check API endpoints.

Tests the /health endpoint which provides basic health monitoring
for the Flask application.
"""

import pytest


class TestHealthAPI:
    """Integration tests for health check endpoints."""

    def test_health_check_returns_ok(self, client):
        """
        GIVEN the Flask application is running
        WHEN GET /health is called
        THEN it returns 200 OK with health status
        """
        # ACT
        response = client.get('/health')

        # ASSERT
        assert response.status_code == 200
        data = response.get_json()
        
        assert 'status' in data
        assert data['status'] == 'ok'
        assert 'service' in data
        assert data['service'] == 'pyscf-front-api'
        assert 'version' in data

    def test_health_check_json_format(self, client):
        """
        GIVEN the Flask application is running
        WHEN GET /health is called
        THEN it returns valid JSON with expected fields
        """
        # ACT
        response = client.get('/health')

        # ASSERT
        assert response.status_code == 200
        assert response.content_type == 'application/json'
        
        data = response.get_json()
        assert isinstance(data, dict)
        assert len(data) >= 3  # At least status, service, version

    def test_health_check_version_field(self, client):
        """
        GIVEN the Flask application is running
        WHEN GET /health is called
        THEN the version field is present (even if 'unknown')
        """
        # ACT
        response = client.get('/health')

        # ASSERT
        data = response.get_json()
        assert 'version' in data
        assert isinstance(data['version'], str)
        # Version should be non-empty string (either version number or 'unknown')
        assert len(data['version']) > 0
