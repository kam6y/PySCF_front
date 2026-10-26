"""
Integration tests for PubChem API endpoints.

Tests the PubChem search and validation endpoints, ensuring proper
handling of various search types and error conditions.
"""

import pytest
from services.exceptions import NotFoundError, ServiceError


class TestPubChemSearchAPI:
    """Integration tests for /api/pubchem/search endpoint."""

    def test_search_by_name_success(self, client, mocker):
        """
        GIVEN PubChemService returns valid XYZ data
        WHEN POST /api/pubchem/search is called with a compound name
        THEN 200 OK is returned with XYZ coordinates
        """
        # ARRANGE
        mock_result = {
            "xyz": "O 0.0000 0.0000 0.1173\nH 0.0000 0.7572 -0.4692\nH 0.0000 -0.7572 -0.4692",
            "cid": 962,
            "name": "Water",
        }
        mock_service = mocker.patch("api.pubchem.get_pubchem_service")
        mock_service.return_value.search_compound.return_value = mock_result

        # ACT
        response = client.post(
            "/api/pubchem/search", json={"query": "water", "searchType": "name"}
        )

        # ASSERT
        assert response.status_code == 200
        data = response.json()
        assert data["success"] is True
        assert "data" in data
        assert data["data"]["xyz"] == mock_result["xyz"]
        assert data["data"]["cid"] == 962

        # Verify service was called correctly
        mock_service.return_value.search_compound.assert_called_once_with(
            "water", "name"
        )

    def test_search_defaults_to_name_when_search_type_omitted(self, client, mocker):
        """
        GIVEN searchType is omitted from the request
        WHEN POST /api/pubchem/search is called
        THEN the service is called with the default name search type
        """
        # ARRANGE
        mock_result = {"xyz": "H 0 0 0", "cid": 123}
        mock_service = mocker.patch("api.pubchem.get_pubchem_service")
        mock_service.return_value.search_compound.return_value = mock_result

        # ACT
        response = client.post("/api/pubchem/search", json={"query": "water"})

        # ASSERT
        assert response.status_code == 200
        mock_service.return_value.search_compound.assert_called_once_with(
            "water", "name"
        )

    @pytest.mark.parametrize("search_type", ["name", "cid", "formula"])
    def test_search_different_types(self, client, mocker, search_type):
        """
        GIVEN PubChemService is configured
        WHEN POST /api/pubchem/search is called with different search types
        THEN the service is called with the correct search type
        """
        # ARRANGE
        mock_result = {"xyz": "H 0 0 0", "cid": 123}
        mock_service = mocker.patch("api.pubchem.get_pubchem_service")
        mock_service.return_value.search_compound.return_value = mock_result

        # ACT
        response = client.post(
            "/api/pubchem/search",
            json={"query": "test_query", "searchType": search_type},
        )

        # ASSERT
        assert response.status_code == 200
        mock_service.return_value.search_compound.assert_called_once_with(
            "test_query", search_type
        )

    def test_search_not_found(self, client, mocker):
        """
        GIVEN PubChemService raises NotFoundError
        WHEN POST /api/pubchem/search is called
        THEN 404 Not Found is returned
        """
        # ARRANGE
        mock_service = mocker.patch("api.pubchem.get_pubchem_service")
        mock_service.return_value.search_compound.side_effect = NotFoundError(
            "Compound not found"
        )

        # ACT
        response = client.post(
            "/api/pubchem/search",
            json={"query": "nonexistent_compound_xyz123", "searchType": "name"},
        )

        # ASSERT
        assert response.status_code == 404
        data = response.json()
        assert data["success"] is False
        assert "error" in data
        assert "not found" in data["error"].lower()

    def test_search_service_error(self, client, mocker):
        """
        GIVEN PubChemService raises ServiceError
        WHEN POST /api/pubchem/search is called
        THEN appropriate error status is returned
        """
        # ARRANGE
        mock_service = mocker.patch("api.pubchem.get_pubchem_service")
        mock_service.return_value.search_compound.side_effect = ServiceError(
            "PubChem API unavailable", status_code=503
        )

        # ACT
        response = client.post(
            "/api/pubchem/search", json={"query": "water", "searchType": "name"}
        )

        # ASSERT
        assert response.status_code == 503
        data = response.json()
        assert data["success"] is False
        # 503 is a curated operational message -- not redacted
        assert data["error"] == "PubChem API unavailable"

    def test_search_missing_required_fields(self, client):
        """
        GIVEN request is missing required fields
        WHEN POST /api/pubchem/search is called
        THEN 400 Bad Request is returned
        """
        # ACT - Missing query (required field)
        response = client.post("/api/pubchem/search", json={"searchType": "name"})

        # ASSERT
        assert response.status_code == 400

    def test_search_invalid_json(self, client):
        """
        GIVEN invalid JSON payload
        WHEN POST /api/pubchem/search is called
        THEN 400 Bad Request is returned
        """
        # ACT
        response = client.post(
            "/api/pubchem/search",
            content=b"invalid json",
            headers={"Content-Type": "application/json"},
        )

        # ASSERT
        assert response.status_code == 400


class TestPubChemValidateAPI:
    """Integration tests for /api/pubchem/validate endpoint."""

    def test_validate_valid_xyz(self, client, mocker):
        """
        GIVEN PubChemService returns validation success
        WHEN POST /api/pubchem/validate is called with valid XYZ
        THEN 200 OK is returned with validation result
        """
        # ARRANGE
        valid_xyz = "H 0 0 0\nH 0 0 0.74"
        mock_result = {"valid": True, "atom_count": 2, "message": "XYZ format is valid"}
        mock_service = mocker.patch("api.pubchem.get_pubchem_service")
        mock_service.return_value.validate_xyz.return_value = mock_result

        # ACT
        response = client.post("/api/pubchem/validate", json={"xyz": valid_xyz})

        # ASSERT
        assert response.status_code == 200
        data = response.json()
        assert data["success"] is True
        assert data["data"]["valid"] is True
        assert data["data"]["atom_count"] == 2

        mock_service.return_value.validate_xyz.assert_called_once_with(valid_xyz)

    def test_validate_invalid_xyz(self, client, mocker):
        """
        GIVEN PubChemService returns validation failure
        WHEN POST /api/pubchem/validate is called with invalid XYZ
        THEN 200 OK is returned with validation error details
        """
        # ARRANGE
        invalid_xyz = "invalid xyz format"
        mock_result = {"valid": False, "message": "Invalid XYZ format"}
        mock_service = mocker.patch("api.pubchem.get_pubchem_service")
        mock_service.return_value.validate_xyz.return_value = mock_result

        # ACT
        response = client.post("/api/pubchem/validate", json={"xyz": invalid_xyz})

        # ASSERT
        assert response.status_code == 200
        data = response.json()
        assert data["success"] is True
        assert data["data"]["valid"] is False

    def test_validate_empty_xyz(self, client):
        """
        GIVEN empty XYZ string
        WHEN POST /api/pubchem/validate is called
        THEN 400 Bad Request is returned
        """
        # ACT
        response = client.post("/api/pubchem/validate", json={"xyz": ""})

        # ASSERT
        assert response.status_code == 400

    def test_validate_missing_xyz_field(self, client):
        """
        GIVEN request is missing xyz field
        WHEN POST /api/pubchem/validate is called
        THEN 400 Bad Request is returned
        """
        # ACT
        response = client.post("/api/pubchem/validate", json={})

        # ASSERT
        assert response.status_code == 400


class TestPubChemExtraFieldRejection:
    """Security tests: extra='forbid' rejects unknown fields at the Pydantic layer."""

    def test_search_rejects_unknown_fields(self, client):
        """
        GIVEN a valid PubChem search payload with an extra unknown field
        WHEN POST /api/pubchem/search is called
        THEN 400 is returned because PubChemSearchRequest has extra='forbid'
        """
        payload = {"query": "water", "searchType": "name", "unknown_field": "x"}

        response = client.post("/api/pubchem/search", json=payload)

        assert response.status_code == 400
        body = response.json()
        assert body["success"] is False
        assert "unknown_field" in body["error"]

    def test_xyz_validate_rejects_unknown_fields(self, client):
        """
        GIVEN a valid XYZ validation payload with an extra unknown field
        WHEN POST /api/pubchem/validate is called
        THEN 400 is returned because XYZValidateRequest has extra='forbid'
        """
        payload = {"xyz": "H 0 0 0\nH 0 0 0.74", "unknown_field": "x"}

        response = client.post("/api/pubchem/validate", json=payload)

        assert response.status_code == 400
        body = response.json()
        assert body["success"] is False
        assert "unknown_field" in body["error"]


class TestPubChemInputLengthLimits:
    """Security tests for input length validation on PubChem endpoints."""

    def test_search_query_over_max_length_returns_400(self, client):
        """
        GIVEN a search query exceeding MAX_QUERY_LENGTH (500)
        WHEN POST /api/pubchem/search is called
        THEN 400 Bad Request is returned
        """
        # ARRANGE
        overlength_query = "a" * 501

        # ACT
        response = client.post(
            "/api/pubchem/search",
            json={
                "query": overlength_query,
                "searchType": "name",
            },
        )

        # ASSERT
        assert response.status_code == 400
        data = response.json()
        assert data["success"] is False

    def test_search_query_at_max_length_is_not_rejected_for_length(
        self, client, mocker
    ):
        """
        GIVEN a search query exactly at MAX_QUERY_LENGTH (500)
        WHEN POST /api/pubchem/search is called
        THEN the request is not rejected for length (200 OK with mocked service)
        """
        # ARRANGE
        boundary_query = "a" * 500
        mock_result = {"xyz": "H 0 0 0", "cid": 123}
        mock_service = mocker.patch("api.pubchem.get_pubchem_service")
        mock_service.return_value.search_compound.return_value = mock_result

        # ACT
        response = client.post(
            "/api/pubchem/search",
            json={
                "query": boundary_query,
                "searchType": "name",
            },
        )

        # ASSERT
        assert response.status_code == 200
        data = response.json()
        assert data["success"] is True

    def test_validate_xyz_over_max_length_returns_400(self, client):
        """
        GIVEN an XYZ string exceeding MAX_XYZ_LENGTH (1_000_000)
        WHEN POST /api/pubchem/validate is called
        THEN 400 Bad Request is returned
        """
        # ARRANGE
        overlength_xyz = "H 0 0 0\n" * 125_001  # > 1_000_000 chars

        # ACT
        response = client.post(
            "/api/pubchem/validate",
            json={
                "xyz": overlength_xyz,
            },
        )

        # ASSERT
        assert response.status_code == 400
        data = response.json()
        assert data["success"] is False
