"""
Integration tests for SMILES Conversion API endpoints.

Tests the SMILES-to-XYZ conversion endpoint, ensuring proper
handling of various SMILES strings and error conditions.
"""

from services.exceptions import ValidationError, ServiceError


class TestSMILESConvertAPI:
    """Integration tests for /api/smiles/convert endpoint."""

    def test_convert_simple_smiles_success(self, client, mocker):
        """
        GIVEN SMILESService returns valid XYZ data
        WHEN POST /api/smiles/convert is called with a valid SMILES string
        THEN 200 OK is returned with XYZ coordinates
        """
        # ARRANGE
        mock_result = {
            "xyz": "O 0.0000 0.0000 0.1173\nH 0.0000 0.7572 -0.4692\nH 0.0000 -0.7572 -0.4692",
            "smiles": "O",
        }
        mock_service = mocker.patch("api.smiles.get_smiles_service")
        mock_service.return_value.convert_smiles.return_value = mock_result

        # ACT
        response = client.post("/api/smiles/convert", json={"smiles": "O"})

        # ASSERT
        assert response.status_code == 200
        data = response.json()
        assert data["success"] is True
        assert "data" in data
        assert "xyz" in data["data"]
        assert "O" in data["data"]["xyz"]
        assert "H" in data["data"]["xyz"]

        # Verify service was called correctly
        mock_service.return_value.convert_smiles.assert_called_once_with("O")

    def test_convert_passes_request_smiles_to_service(self, client, mocker):
        """
        GIVEN SMILES string has leading/trailing whitespace
        WHEN POST /api/smiles/convert is called
        THEN the API delegates validation and trimming to the service layer
        """
        # ARRANGE
        smiles_with_whitespace = "  CCO  "
        mock_result = {"xyz": "C 0 0 0\nC 1 0 0\nO 2 0 0", "smiles": "CCO"}
        mock_service = mocker.patch("api.smiles.get_smiles_service")
        mock_service.return_value.convert_smiles.return_value = mock_result

        # ACT
        response = client.post(
            "/api/smiles/convert", json={"smiles": smiles_with_whitespace}
        )

        # ASSERT
        assert response.status_code == 200
        mock_service.return_value.convert_smiles.assert_called_once_with(
            smiles_with_whitespace
        )

    def test_convert_invalid_smiles(self, client, mocker):
        """
        GIVEN SMILESService raises ValidationError for invalid SMILES
        WHEN POST /api/smiles/convert is called
        THEN 400 Bad Request is returned
        """
        # ARRANGE
        invalid_smiles = "INVALID_SMILES_XXX"
        mock_service = mocker.patch("api.smiles.get_smiles_service")
        mock_service.return_value.convert_smiles.side_effect = ValidationError(
            "Invalid SMILES string"
        )

        # ACT
        response = client.post("/api/smiles/convert", json={"smiles": invalid_smiles})

        # ASSERT
        assert response.status_code == 400
        data = response.json()
        assert data["success"] is False
        assert "error" in data

    def test_convert_empty_smiles(self, client):
        """
        GIVEN empty SMILES string
        WHEN POST /api/smiles/convert is called
        THEN 400 Bad Request is returned
        """
        # ACT
        response = client.post("/api/smiles/convert", json={"smiles": ""})

        # ASSERT
        assert response.status_code == 400

    def test_convert_missing_smiles_field(self, client):
        """
        GIVEN request is missing smiles field
        WHEN POST /api/smiles/convert is called
        THEN 400 Bad Request is returned
        """
        # ACT
        response = client.post("/api/smiles/convert", json={})

        # ASSERT
        assert response.status_code == 400

    def test_convert_service_error(self, client, mocker):
        """
        GIVEN SMILESService raises ServiceError
        WHEN POST /api/smiles/convert is called
        THEN appropriate error status is returned
        """
        # ARRANGE
        mock_service = mocker.patch("api.smiles.get_smiles_service")
        mock_service.return_value.convert_smiles.side_effect = ServiceError(
            "RDKit library error", status_code=500
        )

        # ACT
        response = client.post("/api/smiles/convert", json={"smiles": "CCO"})

        # ASSERT
        assert response.status_code == 500
        data = response.json()
        assert data["success"] is False
        assert data["error"] == "An internal server error occurred."

    def test_convert_whitespace_only_smiles(self, client):
        """
        GIVEN SMILES string contains only whitespace
        WHEN POST /api/smiles/convert is called
        THEN 400 Bad Request is returned
        """
        # ACT
        response = client.post("/api/smiles/convert", json={"smiles": "   "})

        # ASSERT
        assert response.status_code == 400

    def test_convert_invalid_json(self, client):
        """
        GIVEN invalid JSON payload
        WHEN POST /api/smiles/convert is called
        THEN 400 Bad Request is returned
        """
        # ACT
        response = client.post(
            "/api/smiles/convert",
            content=b"invalid json",
            headers={"Content-Type": "application/json"},
        )

        # ASSERT
        assert response.status_code == 400


class TestSMILESExtraFieldRejection:
    """Security tests: extra='forbid' rejects unknown fields at the Pydantic layer."""

    def test_convert_rejects_unknown_fields(self, client):
        """
        GIVEN a valid SMILES convert payload with an extra unknown field
        WHEN POST /api/smiles/convert is called
        THEN 400 is returned because SMILESConvertRequest has extra='forbid'
        """
        payload = {"smiles": "O", "unknown_field": "x"}

        response = client.post("/api/smiles/convert", json=payload)

        assert response.status_code == 400
        body = response.json()
        assert body["success"] is False
        assert "unknown_field" in body["error"]


class TestSMILESInputLengthLimits:
    """Security tests for input length validation on SMILES endpoints."""

    def test_convert_smiles_over_max_length_returns_400(self, client):
        """
        GIVEN a SMILES string exceeding MAX_SMILES_LENGTH (10_000)
        WHEN POST /api/smiles/convert is called
        THEN 400 Bad Request is returned
        """
        # ARRANGE
        overlength_smiles = "C" * 10_001

        # ACT
        response = client.post(
            "/api/smiles/convert",
            json={
                "smiles": overlength_smiles,
            },
        )

        # ASSERT
        assert response.status_code == 400
        data = response.json()
        assert data["success"] is False

    def test_convert_smiles_at_max_length_is_not_rejected_for_length(
        self, client, mocker
    ):
        """
        GIVEN a SMILES string exactly at MAX_SMILES_LENGTH (10_000)
        WHEN POST /api/smiles/convert is called
        THEN the request is not rejected for length (200 OK with mocked service)
        """
        # ARRANGE
        boundary_smiles = "C" * 10_000
        mock_result = {"xyz": "C 0 0 0", "smiles": boundary_smiles}
        mock_service = mocker.patch("api.smiles.get_smiles_service")
        mock_service.return_value.convert_smiles.return_value = mock_result

        # ACT
        response = client.post(
            "/api/smiles/convert",
            json={
                "smiles": boundary_smiles,
            },
        )

        # ASSERT
        assert response.status_code == 200
        data = response.json()
        assert data["success"] is True
