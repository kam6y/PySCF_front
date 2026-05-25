"""
Integration tests for Health Check API endpoints.

Tests the /health endpoint which provides basic health monitoring
for the FastAPI application.
"""

class TestHealthAPI:
    """Integration tests for health check endpoints."""

    def test_health_check_returns_ok(self, client):
        """
        GIVEN the FastAPI application is running
        WHEN GET /health is called
        THEN it returns 200 OK with health status
        """
        # ACT
        response = client.get('/health')

        # ASSERT
        assert response.status_code == 200
        data = response.json()
        
        assert 'status' in data
        assert data['status'] == 'ok'
        assert 'service' in data
        assert data['service'] == 'pyscf-front-api'
        assert 'version' in data
        assert response.headers['content-type'] == 'application/json'
        assert isinstance(data, dict)
        assert len(data) >= 3  # At least status, service, version
        assert isinstance(data['version'], str)
        assert len(data['version']) > 0
