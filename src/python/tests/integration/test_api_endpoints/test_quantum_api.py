"""
Integration tests for Quantum Chemistry API endpoints.

Tests the quantum calculation endpoints including job submission,
monitoring, results retrieval, and orbital/spectrum analysis.
"""

import pytest
from quantum_calc._calculation_repository import CalculationRepository
from services.exceptions import NotFoundError, ServiceError, ValidationError
from services.quantum_service import QuantumService


class TestSupportedParametersAPI:
    """Integration tests for /api/quantum/supported-parameters endpoint."""

    def test_get_supported_parameters_success(self, client, mocker):
        """
        GIVEN QuantumService returns supported parameters
        WHEN GET /api/quantum/supported-parameters is called
        THEN 200 OK is returned with parameter lists
        """
        # ARRANGE
        mock_params = {
            "basis_sets": ["sto-3g", "6-31g", "cc-pvdz"],
            "functionals": ["b3lyp", "pbe0", "m06-2x"],
            "solvents": ["water", "ethanol", "acetone"],
            "calculation_methods": ["HF", "DFT", "MP2", "CCSD"],
        }
        mock_service = mocker.patch("api.quantum.get_quantum_service")
        mock_service.return_value.get_supported_parameters.return_value = mock_params

        # ACT
        response = client.get("/api/quantum/supported-parameters")

        # ASSERT
        assert response.status_code == 200
        data = response.json()
        assert data["success"] is True
        assert "data" in data
        assert "basis_sets" in data["data"]
        assert "functionals" in data["data"]
        assert len(data["data"]["basis_sets"]) > 0


class TestCalculationSubmissionAPI:
    """Integration tests for POST /api/quantum/calculate endpoint."""

    def test_start_calculation_success(self, client, mocker, valid_dft_params):
        """
        GIVEN QuantumService successfully starts a calculation
        WHEN POST /api/quantum/calculate is called with valid parameters
        THEN 202 Accepted is returned with calculation instance
        """
        # ARRANGE
        mock_calc_instance = {
            "id": "calc-123",
            "name": "Test H2 DFT",
            "status": "pending",
            "createdAt": "2024-01-01T00:00:00",
            "parameters": valid_dft_params,
        }
        mock_service = mocker.patch("api.quantum.get_quantum_service")
        mock_service.return_value.start_calculation.return_value = mock_calc_instance

        # ACT
        response = client.post("/api/quantum/calculate", json=valid_dft_params)

        # ASSERT
        assert response.status_code == 202
        data = response.json()
        assert data["success"] is True
        assert "calculation" in data["data"]
        assert data["data"]["calculation"]["id"] == "calc-123"
        assert data["data"]["calculation"]["status"] in [
            "pending",
            "waiting",
            "running",
        ]

    def test_start_hf_calculation(self, client, mocker, valid_hf_params):
        """
        GIVEN valid Hartree-Fock parameters
        WHEN POST /api/quantum/calculate is called
        THEN calculation is started successfully
        """
        # ARRANGE
        mock_calc_instance = {
            "id": "calc-hf-123",
            "name": "Test H2 HF",
            "status": "pending",
            "createdAt": "2024-01-01T00:00:00",
            "parameters": valid_hf_params,
        }
        mock_service = mocker.patch("api.quantum.get_quantum_service")
        mock_service.return_value.start_calculation.return_value = mock_calc_instance

        # ACT
        response = client.post("/api/quantum/calculate", json=valid_hf_params)

        # ASSERT
        assert response.status_code == 202
        data = response.json()
        assert data["success"] is True

    @pytest.mark.parametrize(
        "invalid_field,invalid_value",
        [
            ("charges", "invalid"),
            ("spin", -1),
            ("xyz", ""),
        ],
    )
    def test_start_calculation_invalid_params(
        self, client, valid_dft_params, invalid_field, invalid_value
    ):
        """
        GIVEN invalid calculation parameters
        WHEN POST /api/quantum/calculate is called
        THEN 400 Bad Request is returned
        """
        # ARRANGE
        invalid_params = {**valid_dft_params, invalid_field: invalid_value}

        # ACT
        response = client.post("/api/quantum/calculate", json=invalid_params)

        # ASSERT
        assert response.status_code in [400, 422]  # Bad Request or Unprocessable Entity

    @pytest.mark.parametrize(
        "invalid_field,invalid_value",
        [
            ("charges", -11),
            ("charges", 11),
            ("spin", 11),
            ("cpu_cores", 33),
            ("memory_mb", 511),
            ("memory_mb", 32769),
        ],
    )
    def test_start_calculation_rejects_openapi_numeric_bounds(
        self,
        client,
        valid_dft_params,
        invalid_field,
        invalid_value,
    ):
        """
        GIVEN calculation parameters outside OpenAPI numeric bounds
        WHEN POST /api/quantum/calculate is called
        THEN 400 Bad Request is returned
        """
        # ARRANGE
        invalid_params = {**valid_dft_params, invalid_field: invalid_value}

        # ACT
        response = client.post("/api/quantum/calculate", json=invalid_params)

        # ASSERT
        assert response.status_code == 400
        data = response.json()
        assert data["success"] is False
        assert invalid_field in data["error"]

    def test_start_calculation_missing_required_fields(self, client):
        """
        GIVEN request is missing required fields
        WHEN POST /api/quantum/calculate is called
        THEN 400 Bad Request is returned
        """
        # ACT - Missing most required fields
        response = client.post("/api/quantum/calculate", json={"name": "Test"})

        # ASSERT
        assert response.status_code == 400

    def test_start_calculation_service_error(self, client, mocker, valid_dft_params):
        """
        GIVEN QuantumService raises ServiceError
        WHEN POST /api/quantum/calculate is called
        THEN appropriate error status is returned
        """
        # ARRANGE
        mock_service = mocker.patch("api.quantum.get_quantum_service")
        mock_service.return_value.start_calculation.side_effect = ServiceError(
            "Failed to create calculation directory", status_code=500
        )

        # ACT
        response = client.post("/api/quantum/calculate", json=valid_dft_params)

        # ASSERT
        assert response.status_code == 500
        data = response.json()
        assert data["success"] is False
        assert data["error"] == "An internal server error occurred."

    def test_dft_rejects_casci_parameters(self, client, sample_h2_xyz):
        """
        GIVEN DFT calculation with CASCI-specific parameters (ncas, nelecas)
        WHEN POST /api/quantum/calculate is called
        THEN 400 Bad Request is returned with parameter applicability error
        """
        # ARRANGE
        invalid_params = {
            "name": "Test DFT with invalid params",
            "xyz": sample_h2_xyz,
            "calculation_method": "DFT",
            "basis_function": "sto-3g",
            "exchange_correlation": "b3lyp",
            "ncas": 4,  # Not applicable to DFT
            "nelecas": 4,  # Not applicable to DFT
        }

        # ACT
        response = client.post("/api/quantum/calculate", json=invalid_params)

        # ASSERT
        assert response.status_code == 400
        data = response.json()
        assert data["success"] is False
        assert "error" in data
        error_message = data["error"].lower()
        assert "ncas" in error_message or "not applicable" in error_message

    def test_tddft_rejects_optimize_geometry_true(self, client, sample_h2_xyz):
        """
        GIVEN TDDFT calculation with optimize_geometry=True (disabled parameter)
        WHEN POST /api/quantum/calculate is called
        THEN 400 Bad Request is returned with Pydantic validation error
        """
        # ARRANGE
        invalid_params = {
            "name": "Test TDDFT with invalid params",
            "xyz": sample_h2_xyz,
            "calculation_method": "TDDFT",
            "basis_function": "sto-3g",
            "exchange_correlation": "b3lyp",
            "tddft_nstates": 10,
            "optimize_geometry": True,  # Disabled for TDDFT
        }

        # ACT
        response = client.post("/api/quantum/calculate", json=invalid_params)

        # ASSERT
        assert response.status_code == 400
        data = response.json()
        # Pydantic validation error is now handled by global errorhandler
        assert "error" in data
        assert "optimize_geometry" in data["error"].lower()


class TestCalculationListAPI:
    """Integration tests for GET /api/quantum/calculations endpoint."""

    def test_list_calculations_empty(self, client, mocker):
        """
        GIVEN no calculations exist
        WHEN GET /api/quantum/calculations is called
        THEN 200 OK is returned with empty list
        """
        # ARRANGE
        mock_result = {"calculations": [], "count": 0}
        mock_service = mocker.patch("api.quantum.get_quantum_service")
        mock_service.return_value.list_calculations.return_value = mock_result

        # ACT
        response = client.get("/api/quantum/calculations")

        # ASSERT
        assert response.status_code == 200
        data = response.json()
        assert data["success"] is True
        assert data["data"]["calculations"] == []
        assert data["data"]["count"] == 0

    def test_list_calculations_with_data(self, client, mocker):
        """
        GIVEN multiple calculations exist
        WHEN GET /api/quantum/calculations is called
        THEN 200 OK is returned with calculation list
        """
        # ARRANGE
        mock_calculations = [
            {
                "id": "calc-1",
                "name": "Test 1",
                "status": "completed",
                "date": "2024-01-01T00:00:00",
            },
            {
                "id": "calc-2",
                "name": "Test 2",
                "status": "running",
                "date": "2024-01-02T00:00:00",
            },
            {
                "id": "calc-3",
                "name": "Test 3",
                "status": "error",
                "date": "2024-01-03T00:00:00",
            },
        ]
        mock_result = {"calculations": mock_calculations, "count": 3}
        mock_service = mocker.patch("api.quantum.get_quantum_service")
        mock_service.return_value.list_calculations.return_value = mock_result

        # ACT
        response = client.get("/api/quantum/calculations")

        # ASSERT
        assert response.status_code == 200
        data = response.json()
        assert data["success"] is True
        assert len(data["data"]["calculations"]) == 3
        assert data["data"]["count"] == 3


class TestCalculationDetailsAPI:
    """Integration tests for GET /api/quantum/calculations/<id> endpoint."""

    def test_get_calculation_details_success(self, client, mocker):
        """
        GIVEN calculation exists
        WHEN GET /api/quantum/calculations/<id> is called
        THEN 200 OK is returned with calculation details
        """
        # ARRANGE
        calc_id = "calc-123"
        mock_calc = {
            "calculation": {
                "id": calc_id,
                "name": "Test Calculation",
                "status": "completed",
                "results": {"energy": -1.06},
            }
        }
        mock_service = mocker.patch("api.quantum.get_quantum_service")
        mock_service.return_value.get_calculation_details.return_value = mock_calc

        # ACT
        response = client.get(f"/api/quantum/calculations/{calc_id}")

        # ASSERT
        assert response.status_code == 200
        data = response.json()
        assert data["success"] is True
        assert data["data"]["calculation"]["id"] == calc_id
        assert "results" in data["data"]["calculation"]

    def test_get_calculation_details_not_found(self, client, mocker):
        """
        GIVEN calculation does not exist
        WHEN GET /api/quantum/calculations/<id> is called
        THEN 404 Not Found is returned
        """
        # ARRANGE
        calc_id = "nonexistent-calc"
        mock_service = mocker.patch("api.quantum.get_quantum_service")
        mock_service.return_value.get_calculation_details.side_effect = NotFoundError(
            f"Calculation {calc_id} not found"
        )

        # ACT
        response = client.get(f"/api/quantum/calculations/{calc_id}")

        # ASSERT
        assert response.status_code == 404
        data = response.json()
        assert data["success"] is False


class TestCalculationUpdateAPI:
    """Integration tests for PUT /api/quantum/calculations/<id> endpoint."""

    def test_update_calculation_name_success(self, client, mocker):
        """
        GIVEN calculation exists
        WHEN PUT /api/quantum/calculations/<id> is called with new name
        THEN 200 OK is returned with updated calculation
        """
        # ARRANGE
        calc_id = "calc-123"
        new_name = "Updated Calculation Name"
        mock_result = {
            "calculation": {"id": calc_id, "name": new_name, "status": "completed"}
        }
        mock_service = mocker.patch("api.quantum.get_quantum_service")
        mock_service.return_value.update_calculation.return_value = mock_result

        # ACT
        response = client.put(
            f"/api/quantum/calculations/{calc_id}", json={"name": new_name}
        )

        # ASSERT
        assert response.status_code == 200
        data = response.json()
        assert data["success"] is True
        assert data["data"]["calculation"]["name"] == new_name

        mock_service.return_value.update_calculation.assert_called_once_with(
            calc_id, new_name
        )

    def test_update_calculation_not_found(self, client, mocker):
        """
        GIVEN calculation does not exist
        WHEN PUT /api/quantum/calculations/<id> is called
        THEN 404 Not Found is returned
        """
        # ARRANGE
        calc_id = "nonexistent-calc"
        mock_service = mocker.patch("api.quantum.get_quantum_service")
        mock_service.return_value.update_calculation.side_effect = NotFoundError(
            f"Calculation {calc_id} not found"
        )

        # ACT
        response = client.put(
            f"/api/quantum/calculations/{calc_id}", json={"name": "New Name"}
        )

        # ASSERT
        assert response.status_code == 404


class TestCalculationPauseAPI:
    """Integration tests for POST /api/quantum/calculations/<id>/pause."""

    def test_pause_stale_running_calculation_returns_400_and_recovers_error(
        self,
        client,
        mocker,
        tmp_path,
    ):
        """
        GIVEN status.json says running but the process manager owns no worker
        WHEN the pause endpoint is called
        THEN the API rejects the pause and persists error status
        """
        base_dir = tmp_path / "calculations"
        base_dir.mkdir()
        service = QuantumService()
        service.repository = CalculationRepository(base_dir=str(base_dir))

        calc_id = "stale-running-calc"
        calc_dir = base_dir / calc_id
        calc_dir.mkdir()
        service.repository.save_calculation_parameters(
            str(calc_dir),
            {"name": "Stale Running Calc", "created_at": "2026-05-20T00:00:00"},
        )
        service.repository.save_calculation_status(str(calc_dir), "running")

        process_manager = mocker.Mock()
        process_manager.get_active_calculations.return_value = []
        process_manager.get_queued_calculations.return_value = []
        process_manager.pause_calculation.return_value = True
        mocker.patch(
            "services.calculation_service_context.get_process_manager",
            return_value=process_manager,
        )
        mocker.patch("api.quantum.get_quantum_service", return_value=service)

        response = client.post(f"/api/quantum/calculations/{calc_id}/pause")

        assert response.status_code == 400
        data = response.json()
        assert data["success"] is False
        assert "status: error" in data["error"]
        process_manager.pause_calculation.assert_not_called()
        assert service.repository.read_calculation_status_details(str(calc_dir)) == (
            "error",
            None,
        )
        assert service.repository.read_calculation_results(str(calc_dir)) == {
            "error": QuantumService.RESTART_INTERRUPTED_MESSAGE,
        }


class TestCalculationDeletionAPI:
    """Integration tests for DELETE /api/quantum/calculations/<id> endpoint."""

    def test_delete_calculation_success(self, client, mocker):
        """
        GIVEN calculation exists
        WHEN DELETE /api/quantum/calculations/<id> is called
        THEN 200 OK is returned with deletion confirmation
        """
        # ARRANGE
        calc_id = "calc-123"
        mock_result = {
            "deleted_id": calc_id,
            "message": f'Calculation "{calc_id}" has been deleted successfully',
        }
        mock_service = mocker.patch("api.quantum.get_quantum_service")
        mock_service.return_value.delete_calculation.return_value = mock_result

        # ACT
        response = client.delete(f"/api/quantum/calculations/{calc_id}")

        # ASSERT
        assert response.status_code == 200
        data = response.json()
        assert data["success"] is True
        assert data["data"]["deleted_id"] == calc_id

    def test_delete_calculation_not_found(self, client, mocker):
        """
        GIVEN calculation does not exist
        WHEN DELETE /api/quantum/calculations/<id> is called
        THEN 404 Not Found is returned
        """
        # ARRANGE
        calc_id = "nonexistent-calc"
        mock_service = mocker.patch("api.quantum.get_quantum_service")
        mock_service.return_value.delete_calculation.side_effect = NotFoundError(
            f"Calculation {calc_id} not found"
        )

        # ACT
        response = client.delete(f"/api/quantum/calculations/{calc_id}")

        # ASSERT
        assert response.status_code == 404

    def test_delete_calculation_validation_error_returns_400(self, client, mocker):
        """
        GIVEN service rejects deletion for a non-terminal calculation
        WHEN DELETE /api/quantum/calculations/<id> is called
        THEN 400 Bad Request is returned
        """
        # ARRANGE
        calc_id = "running-calc"
        mock_service = mocker.patch("api.quantum.get_quantum_service")
        mock_service.return_value.delete_calculation.side_effect = ValidationError(
            f'Cannot delete calculation "{calc_id}" while it is running.'
        )

        # ACT
        response = client.delete(f"/api/quantum/calculations/{calc_id}")

        # ASSERT
        assert response.status_code == 400
        data = response.json()
        assert data["success"] is False
        assert "Cannot delete calculation" in data["error"]

    def test_delete_calculation_rejects_encoded_parent_directory_and_keeps_sentinel(
        self,
        client,
        mocker,
        tmp_path,
    ):
        """
        GIVEN the calculations base has a parent sentinel directory
        WHEN DELETE receives an encoded parent-directory calculation ID
        THEN it returns 400 and does not remove data outside the base directory
        """
        base_dir = tmp_path / "calculations"
        base_dir.mkdir()
        sentinel_dir = tmp_path / "sentinel"
        sentinel_dir.mkdir()
        sentinel_file = sentinel_dir / "keep.txt"
        sentinel_file.write_text("must remain")

        service = QuantumService()
        service.repository = CalculationRepository(base_dir=str(base_dir))
        process_manager = mocker.Mock()
        process_manager.get_active_calculations.return_value = []
        process_manager.get_queued_calculations.return_value = []
        mocker.patch(
            "services.calculation_service_context.get_process_manager",
            return_value=process_manager,
        )
        mocker.patch("api.quantum.get_quantum_service", return_value=service)

        response = client.delete("/api/quantum/calculations/%2e%2e")

        assert response.status_code == 400
        assert base_dir.is_dir()
        assert sentinel_dir.is_dir()
        assert sentinel_file.read_text() == "must remain"


class TestCalculationIdValidationAPI:
    """Integration tests for request-derived calculation ID validation."""

    @pytest.mark.parametrize(
        ("method", "path", "json_body"),
        [
            ("get", "/api/quantum/calculations/a%5Cb", None),
            ("put", "/api/quantum/calculations/a%5Cb", {"name": "New Name"}),
            ("delete", "/api/quantum/calculations/a%5Cb", None),
            ("post", "/api/quantum/calculations/a%5Cb/pause", None),
            ("post", "/api/quantum/calculations/a%5Cb/resume", None),
            ("get", "/api/quantum/calculations/a%5Cb/orbitals", None),
            ("get", "/api/quantum/calculations/a%5Cb/orbitals/1/cube", None),
            ("get", "/api/quantum/calculations/a%5Cb/orbitals/cube-files", None),
            ("delete", "/api/quantum/calculations/a%5Cb/orbitals/cube-files", None),
            ("get", "/api/quantum/calculations/a%5Cb/ir-spectrum", None),
        ],
    )
    def test_calculation_id_endpoints_reject_path_separator_ids(
        self,
        client,
        mocker,
        tmp_path,
        method,
        path,
        json_body,
    ):
        """
        GIVEN a request-derived calculation ID contains a path separator
        WHEN any calculation-ID endpoint is called
        THEN the service validation returns 400
        """
        service = QuantumService()
        service.repository = CalculationRepository(
            base_dir=str(tmp_path / "calculations")
        )
        mocker.patch("api.quantum.get_quantum_service", return_value=service)

        request_method = getattr(client, method)
        kwargs = {"json": json_body} if json_body is not None else {}
        response = request_method(path, **kwargs)

        assert response.status_code == 400
        data = response.json()
        assert data["success"] is False
        assert "Invalid calculation ID" in data["error"]


class TestMolecularOrbitalsAPI:
    """Integration tests for orbital-related endpoints."""

    def test_get_orbitals_success(self, client, mocker):
        """
        GIVEN calculation has orbital data
        WHEN GET /api/quantum/calculations/<id>/orbitals is called
        THEN 200 OK is returned with orbital information
        """
        # ARRANGE
        calc_id = "calc-123"
        mock_orbitals = {
            "homo_index": 4,
            "lumo_index": 5,
            "orbitals": [
                {"index": 0, "energy": -10.5, "occupancy": 2},
                {"index": 1, "energy": -8.2, "occupancy": 2},
            ],
        }
        mock_service = mocker.patch("api.quantum.get_quantum_service")
        mock_service.return_value.get_molecular_orbitals.return_value = mock_orbitals

        # ACT
        response = client.get(f"/api/quantum/calculations/{calc_id}/orbitals")

        # ASSERT
        assert response.status_code == 200
        data = response.json()
        assert data["success"] is True
        assert "homo_index" in data["data"]
        assert "orbitals" in data["data"]

    def test_generate_orbital_cube_success(self, client, mocker):
        """
        GIVEN calculation exists with orbital data
        WHEN GET /api/quantum/calculations/<id>/orbitals/<index>/cube is called
        THEN 200 OK is returned with CUBE file data
        """
        # ARRANGE
        calc_id = "calc-123"
        orbital_index = 5
        mock_cube_data = {
            "orbital_index": orbital_index,
            "cube_file_content": "CUBE file content here...",
            "generation_info": {"grid_size": 80},
        }
        mock_service = mocker.patch("api.quantum.get_quantum_service")
        mock_service.return_value.generate_orbital_cube.return_value = mock_cube_data

        # ACT
        response = client.get(
            f"/api/quantum/calculations/{calc_id}/orbitals/{orbital_index}/cube"
        )

        # ASSERT
        assert response.status_code == 200
        data = response.json()
        assert data["success"] is True
        assert data["data"]["orbital_index"] == orbital_index

    def test_generate_orbital_cube_with_params(self, client, mocker):
        """
        GIVEN custom CUBE generation parameters
        WHEN GET with query parameters is called
        THEN service receives the custom parameters
        """
        # ARRANGE
        calc_id = "calc-123"
        orbital_index = 5
        mock_cube_data = {"orbital_index": orbital_index}
        mock_service = mocker.patch("api.quantum.get_quantum_service")
        mock_service.return_value.generate_orbital_cube.return_value = mock_cube_data

        # ACT
        response = client.get(
            f"/api/quantum/calculations/{calc_id}/orbitals/{orbital_index}/cube"
            f"?gridSize=100&isovaluePos=0.05&isovalueNeg=-0.05"
        )

        # ASSERT
        assert response.status_code == 200
        mock_service.return_value.generate_orbital_cube.assert_called_once_with(
            calc_id, orbital_index, grid_size=100, isovalue_pos=0.05, isovalue_neg=-0.05
        )

    def test_generate_orbital_cube_validation_error_returns_400(self, client, mocker):
        """
        GIVEN service rejects out-of-range CUBE query parameters
        WHEN GET /orbitals/<index>/cube is called
        THEN 400 Bad Request is returned
        """
        # ARRANGE
        calc_id = "calc-123"
        orbital_index = 5
        mock_service = mocker.patch("api.quantum.get_quantum_service")
        mock_service.return_value.generate_orbital_cube.side_effect = ValidationError(
            "grid_size must be between 40 and 120."
        )

        # ACT
        response = client.get(
            f"/api/quantum/calculations/{calc_id}/orbitals/{orbital_index}/cube?gridSize=121"
        )

        # ASSERT
        assert response.status_code == 400
        data = response.json()
        assert data["success"] is False
        assert "grid_size" in data["error"]

    def test_generate_orbital_cube_invalid_orbital_index_returns_400(
        self, client, mocker
    ):
        """
        GIVEN service rejects an unavailable orbital index
        WHEN GET /orbitals/<index>/cube is called
        THEN 400 Bad Request is returned
        """
        # ARRANGE
        calc_id = "calc-123"
        orbital_index = 999
        mock_service = mocker.patch("api.quantum.get_quantum_service")
        mock_service.return_value.generate_orbital_cube.side_effect = ValidationError(
            "Invalid orbital index: 999. Available range: 0-5"
        )

        # ACT
        response = client.get(
            f"/api/quantum/calculations/{calc_id}/orbitals/{orbital_index}/cube"
        )

        # ASSERT
        assert response.status_code == 400
        data = response.json()
        assert data["success"] is False
        assert "Invalid orbital index" in data["error"]

    def test_list_cube_files_success(self, client, mocker):
        """
        GIVEN calculation has CUBE files
        WHEN GET /api/quantum/calculations/<id>/orbitals/cube-files is called
        THEN 200 OK is returned with file list
        """
        # ARRANGE
        calc_id = "calc-123"
        mock_files = {
            "cube_files": [
                {"orbital_index": 4, "filename": "orbital_4.cube"},
                {"orbital_index": 5, "filename": "orbital_5.cube"},
            ],
            "total_files": 2,
        }
        mock_service = mocker.patch("api.quantum.get_quantum_service")
        mock_service.return_value.list_cube_files.return_value = mock_files

        # ACT
        response = client.get(
            f"/api/quantum/calculations/{calc_id}/orbitals/cube-files"
        )

        # ASSERT
        assert response.status_code == 200
        data = response.json()
        assert data["success"] is True
        assert data["data"]["total_files"] == 2

    def test_delete_cube_files_all(self, client, mocker):
        """
        GIVEN calculation has CUBE files
        WHEN DELETE /api/quantum/calculations/<id>/orbitals/cube-files is called
        THEN all CUBE files are deleted
        """
        # ARRANGE
        calc_id = "calc-123"
        mock_result = {"deleted_files": 5, "message": "All CUBE files deleted"}
        mock_service = mocker.patch("api.quantum.get_quantum_service")
        mock_service.return_value.delete_cube_files.return_value = mock_result

        # ACT
        response = client.delete(
            f"/api/quantum/calculations/{calc_id}/orbitals/cube-files"
        )

        # ASSERT
        assert response.status_code == 200
        data = response.json()
        assert data["success"] is True
        assert data["data"]["deleted_files"] == 5

    def test_delete_cube_files_specific_orbital(self, client, mocker):
        """
        GIVEN orbital_index query parameter
        WHEN DELETE is called
        THEN only that orbital's CUBE file is deleted
        """
        # ARRANGE
        calc_id = "calc-123"
        orbital_index = 5
        mock_result = {"deleted_files": 1, "orbital_index": orbital_index}
        mock_service = mocker.patch("api.quantum.get_quantum_service")
        mock_service.return_value.delete_cube_files.return_value = mock_result

        # ACT
        response = client.delete(
            f"/api/quantum/calculations/{calc_id}/orbitals/cube-files"
            f"?orbital_index={orbital_index}"
        )

        # ASSERT
        assert response.status_code == 200
        mock_service.return_value.delete_cube_files.assert_called_once_with(
            calc_id, orbital_index
        )


class TestIRSpectrumAPI:
    """Integration tests for IR spectrum generation endpoint."""

    def test_generate_ir_spectrum_success(self, client, mocker):
        """
        GIVEN calculation has frequency data
        WHEN GET /api/quantum/calculations/<id>/ir-spectrum is called
        THEN 200 OK is returned with spectrum data
        """
        # ARRANGE
        calc_id = "calc-123"
        mock_spectrum = {
            "spectrum": {"x": [400, 500, 600], "y": [0.1, 0.5, 0.2]},
            "plot_image_base64": "base64encodedimage...",
        }
        mock_service = mocker.patch("api.quantum.get_quantum_service")
        mock_service.return_value.generate_ir_spectrum.return_value = mock_spectrum

        # ACT
        response = client.get(f"/api/quantum/calculations/{calc_id}/ir-spectrum")

        # ASSERT
        assert response.status_code == 200
        data = response.json()
        assert data["success"] is True
        assert "spectrum" in data["data"]

    def test_generate_ir_spectrum_with_custom_params(self, client, mocker):
        """
        GIVEN custom spectrum generation parameters
        WHEN GET with query parameters is called
        THEN service receives the custom parameters
        """
        # ARRANGE
        calc_id = "calc-123"
        mock_spectrum = {"spectrum": {}}
        mock_service = mocker.patch("api.quantum.get_quantum_service")
        mock_service.return_value.generate_ir_spectrum.return_value = mock_spectrum

        # ACT
        response = client.get(
            f"/api/quantum/calculations/{calc_id}/ir-spectrum"
            f"?broadening_fwhm=50&x_min=500&x_max=3500"
        )

        # ASSERT
        assert response.status_code == 200
        mock_service.return_value.generate_ir_spectrum.assert_called_once_with(
            calc_id,
            broadening_fwhm=50.0,
            x_min=500.0,
            x_max=3500.0,
            show_peaks=True,  # Default value when not specified
        )

    def test_generate_ir_spectrum_not_found(self, client, mocker):
        """
        GIVEN calculation does not have frequency data
        WHEN GET /api/quantum/calculations/<id>/ir-spectrum is called
        THEN 404 Not Found is returned
        """
        # ARRANGE
        calc_id = "calc-no-freq"
        mock_service = mocker.patch("api.quantum.get_quantum_service")
        mock_service.return_value.generate_ir_spectrum.side_effect = NotFoundError(
            "No frequency data available"
        )

        # ACT
        response = client.get(f"/api/quantum/calculations/{calc_id}/ir-spectrum")

        # ASSERT
        assert response.status_code == 404


class TestCalculationStatusAPI:
    """Integration tests for GET /api/quantum/status endpoint."""

    def test_get_calculation_system_status(self, client, mocker):
        """
        GIVEN quantum calculation system is running
        WHEN GET /api/quantum/status is called
        THEN 200 OK is returned with system status
        """
        # ARRANGE
        mock_status = {
            "process_pool": {"active": True, "workers": 4},
            "system": {"cpu_count": 8, "memory_available": 16000},
        }
        mock_service = mocker.patch("api.quantum.get_quantum_service")
        mock_service.return_value.get_calculation_status.return_value = mock_status

        # ACT
        response = client.get("/api/quantum/status")

        # ASSERT
        assert response.status_code == 200
        data = response.json()
        assert data["success"] is True
        assert "process_pool" in data["data"]


class TestQuantumInputLengthLimits:
    """Security tests for input length validation on quantum calculation endpoints."""

    @pytest.mark.parametrize(
        "field_name",
        [
            "basis_function",
            "solvent",
            "exchange_correlation",
            "auxiliary_basis",
        ],
    )
    def test_short_field_over_max_length_returns_400(
        self, client, valid_dft_params, field_name
    ):
        """
        GIVEN a short identifier field exceeding MAX_SHORT_FIELD_LENGTH (200)
        WHEN POST /api/quantum/calculate is called
        THEN 400 Bad Request is returned with the field name in the error message
        """
        # ARRANGE
        overlength_params = {**valid_dft_params, field_name: "x" * 201}

        # ACT
        response = client.post("/api/quantum/calculate", json=overlength_params)

        # ASSERT
        assert response.status_code == 400
        data = response.json()
        assert data["success"] is False
        assert field_name in data["error"]

    def test_short_field_at_max_length_is_not_rejected_for_length(
        self, client, mocker, valid_dft_params
    ):
        """
        GIVEN basis_function exactly at MAX_SHORT_FIELD_LENGTH (200)
        WHEN POST /api/quantum/calculate is called
        THEN the request is not rejected for length (202 Accepted with mocked service)
        """
        # ARRANGE
        boundary_params = {**valid_dft_params, "basis_function": "x" * 200}
        mock_calc_instance = {
            "id": "calc-boundary",
            "name": "Boundary Test",
            "status": "pending",
            "createdAt": "2024-01-01T00:00:00",
            "parameters": boundary_params,
        }
        mock_service = mocker.patch("api.quantum.get_quantum_service")
        mock_service.return_value.start_calculation.return_value = mock_calc_instance

        # ACT
        response = client.post("/api/quantum/calculate", json=boundary_params)

        # ASSERT
        assert response.status_code == 202
        data = response.json()
        assert data["success"] is True

    def test_xyz_over_max_length_returns_400(self, client, valid_dft_params):
        """
        GIVEN xyz data exceeding MAX_XYZ_LENGTH (1_000_000)
        WHEN POST /api/quantum/calculate is called
        THEN 400 Bad Request is returned with xyz mentioned in the error message
        """
        # ARRANGE
        overlength_params = {**valid_dft_params, "xyz": "H 0 0 0\n" * 125_001}

        # ACT
        response = client.post("/api/quantum/calculate", json=overlength_params)

        # ASSERT
        assert response.status_code == 400
        data = response.json()
        assert data["success"] is False
        assert "xyz" in data["error"].lower()
