"""
Integration tests for Quantum Chemistry API endpoints.

Tests the quantum calculation endpoints including job submission,
monitoring, results retrieval, and orbital/spectrum analysis.
"""

import pytest
from services.exceptions import NotFoundError, ValidationError, ServiceError


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
            'basis_sets': ['sto-3g', '6-31g', 'cc-pvdz'],
            'functionals': ['b3lyp', 'pbe0', 'm06-2x'],
            'solvents': ['water', 'ethanol', 'acetone'],
            'calculation_methods': ['HF', 'DFT', 'MP2', 'CCSD']
        }
        mock_service = mocker.patch('api.quantum.get_quantum_service')
        mock_service.return_value.get_supported_parameters.return_value = mock_params

        # ACT
        response = client.get('/api/quantum/supported-parameters')

        # ASSERT
        assert response.status_code == 200
        data = response.get_json()
        assert data['success'] is True
        assert 'data' in data
        assert 'basis_sets' in data['data']
        assert 'functionals' in data['data']
        assert len(data['data']['basis_sets']) > 0


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
            'id': 'calc-123',
            'name': 'Test H2 DFT',
            'status': 'pending',
            'createdAt': '2024-01-01T00:00:00',
            'parameters': valid_dft_params
        }
        mock_service = mocker.patch('api.quantum.get_quantum_service')
        mock_service.return_value.start_calculation.return_value = mock_calc_instance

        # ACT
        response = client.post('/api/quantum/calculate', json=valid_dft_params)

        # ASSERT
        assert response.status_code == 202
        data = response.get_json()
        assert data['success'] is True
        assert 'calculation' in data['data']
        assert data['data']['calculation']['id'] == 'calc-123'
        assert data['data']['calculation']['status'] in ['pending', 'waiting', 'running']

    def test_start_hf_calculation(self, client, mocker, valid_hf_params):
        """
        GIVEN valid Hartree-Fock parameters
        WHEN POST /api/quantum/calculate is called
        THEN calculation is started successfully
        """
        # ARRANGE
        mock_calc_instance = {
            'id': 'calc-hf-123',
            'name': 'Test H2 HF',
            'status': 'pending',
            'createdAt': '2024-01-01T00:00:00',
            'parameters': valid_hf_params
        }
        mock_service = mocker.patch('api.quantum.get_quantum_service')
        mock_service.return_value.start_calculation.return_value = mock_calc_instance

        # ACT
        response = client.post('/api/quantum/calculate', json=valid_hf_params)

        # ASSERT
        assert response.status_code == 202
        data = response.get_json()
        assert data['success'] is True

    @pytest.mark.parametrize("invalid_field,invalid_value", [
        ('charges', 'invalid'),
        ('spin', -1),
        ('xyz', ''),
    ])
    def test_start_calculation_invalid_params(self, client, mocker, valid_dft_params, invalid_field, invalid_value):
        """
        GIVEN invalid calculation parameters
        WHEN POST /api/quantum/calculate is called
        THEN 400 Bad Request is returned
        """
        # ARRANGE
        invalid_params = {**valid_dft_params, invalid_field: invalid_value}

        # ACT
        response = client.post('/api/quantum/calculate', json=invalid_params)

        # ASSERT
        assert response.status_code in [400, 422]  # Bad Request or Unprocessable Entity

    def test_start_calculation_missing_required_fields(self, client):
        """
        GIVEN request is missing required fields
        WHEN POST /api/quantum/calculate is called
        THEN 400 Bad Request is returned
        """
        # ACT - Missing most required fields
        response = client.post('/api/quantum/calculate', json={
            'name': 'Test'
        })

        # ASSERT
        assert response.status_code == 400

    def test_start_calculation_service_error(self, client, mocker, valid_dft_params):
        """
        GIVEN QuantumService raises ServiceError
        WHEN POST /api/quantum/calculate is called
        THEN appropriate error status is returned
        """
        # ARRANGE
        mock_service = mocker.patch('api.quantum.get_quantum_service')
        mock_service.return_value.start_calculation.side_effect = ServiceError(
            "Failed to create calculation directory",
            status_code=500
        )

        # ACT
        response = client.post('/api/quantum/calculate', json=valid_dft_params)

        # ASSERT
        assert response.status_code == 500
        data = response.get_json()
        assert data['success'] is False


class TestCalculationListAPI:
    """Integration tests for GET /api/quantum/calculations endpoint."""

    def test_list_calculations_empty(self, client, mocker):
        """
        GIVEN no calculations exist
        WHEN GET /api/quantum/calculations is called
        THEN 200 OK is returned with empty list
        """
        # ARRANGE
        mock_result = {'calculations': [], 'count': 0}
        mock_service = mocker.patch('api.quantum.get_quantum_service')
        mock_service.return_value.list_calculations.return_value = mock_result

        # ACT
        response = client.get('/api/quantum/calculations')

        # ASSERT
        assert response.status_code == 200
        data = response.get_json()
        assert data['success'] is True
        assert data['data']['calculations'] == []
        assert data['data']['count'] == 0

    def test_list_calculations_with_data(self, client, mocker):
        """
        GIVEN multiple calculations exist
        WHEN GET /api/quantum/calculations is called
        THEN 200 OK is returned with calculation list
        """
        # ARRANGE
        mock_calculations = [
            {'id': 'calc-1', 'name': 'Test 1', 'status': 'completed'},
            {'id': 'calc-2', 'name': 'Test 2', 'status': 'running'},
            {'id': 'calc-3', 'name': 'Test 3', 'status': 'error'}
        ]
        mock_result = {'calculations': mock_calculations, 'count': 3}
        mock_service = mocker.patch('api.quantum.get_quantum_service')
        mock_service.return_value.list_calculations.return_value = mock_result

        # ACT
        response = client.get('/api/quantum/calculations')

        # ASSERT
        assert response.status_code == 200
        data = response.get_json()
        assert data['success'] is True
        assert len(data['data']['calculations']) == 3
        assert data['data']['count'] == 3


class TestCalculationDetailsAPI:
    """Integration tests for GET /api/quantum/calculations/<id> endpoint."""

    def test_get_calculation_details_success(self, client, mocker):
        """
        GIVEN calculation exists
        WHEN GET /api/quantum/calculations/<id> is called
        THEN 200 OK is returned with calculation details
        """
        # ARRANGE
        calc_id = 'calc-123'
        mock_calc = {
            'calculation': {
                'id': calc_id,
                'name': 'Test Calculation',
                'status': 'completed',
                'results': {'energy': -1.06}
            }
        }
        mock_service = mocker.patch('api.quantum.get_quantum_service')
        mock_service.return_value.get_calculation_details.return_value = mock_calc

        # ACT
        response = client.get(f'/api/quantum/calculations/{calc_id}')

        # ASSERT
        assert response.status_code == 200
        data = response.get_json()
        assert data['success'] is True
        assert data['data']['calculation']['id'] == calc_id
        assert 'results' in data['data']['calculation']

    def test_get_calculation_details_not_found(self, client, mocker):
        """
        GIVEN calculation does not exist
        WHEN GET /api/quantum/calculations/<id> is called
        THEN 404 Not Found is returned
        """
        # ARRANGE
        calc_id = 'nonexistent-calc'
        mock_service = mocker.patch('api.quantum.get_quantum_service')
        mock_service.return_value.get_calculation_details.side_effect = NotFoundError(f"Calculation {calc_id} not found")

        # ACT
        response = client.get(f'/api/quantum/calculations/{calc_id}')

        # ASSERT
        assert response.status_code == 404
        data = response.get_json()
        assert data['success'] is False


class TestCalculationUpdateAPI:
    """Integration tests for PUT /api/quantum/calculations/<id> endpoint."""

    def test_update_calculation_name_success(self, client, mocker):
        """
        GIVEN calculation exists
        WHEN PUT /api/quantum/calculations/<id> is called with new name
        THEN 200 OK is returned with updated calculation
        """
        # ARRANGE
        calc_id = 'calc-123'
        new_name = 'Updated Calculation Name'
        mock_result = {
            'calculation': {
                'id': calc_id,
                'name': new_name,
                'status': 'completed'
            }
        }
        mock_service = mocker.patch('api.quantum.get_quantum_service')
        mock_service.return_value.update_calculation.return_value = mock_result

        # ACT
        response = client.put(f'/api/quantum/calculations/{calc_id}', json={
            'name': new_name
        })

        # ASSERT
        assert response.status_code == 200
        data = response.get_json()
        assert data['success'] is True
        assert data['data']['calculation']['name'] == new_name
        
        mock_service.return_value.update_calculation.assert_called_once_with(calc_id, new_name)

    def test_update_calculation_not_found(self, client, mocker):
        """
        GIVEN calculation does not exist
        WHEN PUT /api/quantum/calculations/<id> is called
        THEN 404 Not Found is returned
        """
        # ARRANGE
        calc_id = 'nonexistent-calc'
        mock_service = mocker.patch('api.quantum.get_quantum_service')
        mock_service.return_value.update_calculation.side_effect = NotFoundError(f"Calculation {calc_id} not found")

        # ACT
        response = client.put(f'/api/quantum/calculations/{calc_id}', json={
            'name': 'New Name'
        })

        # ASSERT
        assert response.status_code == 404


class TestCalculationDeletionAPI:
    """Integration tests for DELETE /api/quantum/calculations/<id> endpoint."""

    def test_delete_calculation_success(self, client, mocker):
        """
        GIVEN calculation exists
        WHEN DELETE /api/quantum/calculations/<id> is called
        THEN 200 OK is returned with deletion confirmation
        """
        # ARRANGE
        calc_id = 'calc-123'
        mock_result = {
            'deleted_id': calc_id,
            'message': f'Calculation "{calc_id}" has been deleted successfully'
        }
        mock_service = mocker.patch('api.quantum.get_quantum_service')
        mock_service.return_value.delete_calculation.return_value = mock_result

        # ACT
        response = client.delete(f'/api/quantum/calculations/{calc_id}')

        # ASSERT
        assert response.status_code == 200
        data = response.get_json()
        assert data['success'] is True
        assert data['data']['deleted_id'] == calc_id

    def test_delete_calculation_not_found(self, client, mocker):
        """
        GIVEN calculation does not exist
        WHEN DELETE /api/quantum/calculations/<id> is called
        THEN 404 Not Found is returned
        """
        # ARRANGE
        calc_id = 'nonexistent-calc'
        mock_service = mocker.patch('api.quantum.get_quantum_service')
        mock_service.return_value.delete_calculation.side_effect = NotFoundError(f"Calculation {calc_id} not found")

        # ACT
        response = client.delete(f'/api/quantum/calculations/{calc_id}')

        # ASSERT
        assert response.status_code == 404


class TestCalculationCancellationAPI:
    """Integration tests for POST /api/quantum/calculations/<id>/cancel endpoint."""

    def test_cancel_calculation_success(self, client, mocker):
        """
        GIVEN calculation is running
        WHEN POST /api/quantum/calculations/<id>/cancel is called
        THEN 200 OK is returned with cancellation confirmation
        """
        # ARRANGE
        calc_id = 'calc-123'
        mock_result = {
            'calculation_id': calc_id,
            'message': f'Calculation "{calc_id}" has been cancelled successfully'
        }
        mock_service = mocker.patch('api.quantum.get_quantum_service')
        mock_service.return_value.cancel_calculation.return_value = mock_result

        # ACT
        response = client.post(f'/api/quantum/calculations/{calc_id}/cancel')

        # ASSERT
        assert response.status_code == 200
        data = response.get_json()
        assert data['success'] is True

    def test_cancel_calculation_not_found(self, client, mocker):
        """
        GIVEN calculation does not exist
        WHEN POST /api/quantum/calculations/<id>/cancel is called
        THEN 404 Not Found is returned
        """
        # ARRANGE
        calc_id = 'nonexistent-calc'
        mock_service = mocker.patch('api.quantum.get_quantum_service')
        mock_service.return_value.cancel_calculation.side_effect = NotFoundError(f"Calculation {calc_id} not found")

        # ACT
        response = client.post(f'/api/quantum/calculations/{calc_id}/cancel')

        # ASSERT
        assert response.status_code == 404


class TestMolecularOrbitalsAPI:
    """Integration tests for orbital-related endpoints."""

    def test_get_orbitals_success(self, client, mocker):
        """
        GIVEN calculation has orbital data
        WHEN GET /api/quantum/calculations/<id>/orbitals is called
        THEN 200 OK is returned with orbital information
        """
        # ARRANGE
        calc_id = 'calc-123'
        mock_orbitals = {
            'homo_index': 4,
            'lumo_index': 5,
            'orbitals': [
                {'index': 0, 'energy': -10.5, 'occupancy': 2},
                {'index': 1, 'energy': -8.2, 'occupancy': 2},
            ]
        }
        mock_service = mocker.patch('api.quantum.get_quantum_service')
        mock_service.return_value.get_molecular_orbitals.return_value = mock_orbitals

        # ACT
        response = client.get(f'/api/quantum/calculations/{calc_id}/orbitals')

        # ASSERT
        assert response.status_code == 200
        data = response.get_json()
        assert data['success'] is True
        assert 'homo_index' in data['data']
        assert 'orbitals' in data['data']

    def test_generate_orbital_cube_success(self, client, mocker):
        """
        GIVEN calculation exists with orbital data
        WHEN GET /api/quantum/calculations/<id>/orbitals/<index>/cube is called
        THEN 200 OK is returned with CUBE file data
        """
        # ARRANGE
        calc_id = 'calc-123'
        orbital_index = 5
        mock_cube_data = {
            'orbital_index': orbital_index,
            'cube_file_content': 'CUBE file content here...',
            'generation_info': {'grid_size': 80}
        }
        mock_service = mocker.patch('api.quantum.get_quantum_service')
        mock_service.return_value.generate_orbital_cube.return_value = mock_cube_data

        # ACT
        response = client.get(f'/api/quantum/calculations/{calc_id}/orbitals/{orbital_index}/cube')

        # ASSERT
        assert response.status_code == 200
        data = response.get_json()
        assert data['success'] is True
        assert data['data']['orbital_index'] == orbital_index

    def test_generate_orbital_cube_with_params(self, client, mocker):
        """
        GIVEN custom CUBE generation parameters
        WHEN GET with query parameters is called
        THEN service receives the custom parameters
        """
        # ARRANGE
        calc_id = 'calc-123'
        orbital_index = 5
        mock_cube_data = {'orbital_index': orbital_index}
        mock_service = mocker.patch('api.quantum.get_quantum_service')
        mock_service.return_value.generate_orbital_cube.return_value = mock_cube_data

        # ACT
        response = client.get(
            f'/api/quantum/calculations/{calc_id}/orbitals/{orbital_index}/cube'
            f'?gridSize=100&isovaluePos=0.05&isovalueNeg=-0.05'
        )

        # ASSERT
        assert response.status_code == 200
        mock_service.return_value.generate_orbital_cube.assert_called_once_with(
            calc_id,
            orbital_index,
            grid_size=100,
            isovalue_pos=0.05,
            isovalue_neg=-0.05
        )

    def test_list_cube_files_success(self, client, mocker):
        """
        GIVEN calculation has CUBE files
        WHEN GET /api/quantum/calculations/<id>/orbitals/cube-files is called
        THEN 200 OK is returned with file list
        """
        # ARRANGE
        calc_id = 'calc-123'
        mock_files = {
            'cube_files': [
                {'orbital_index': 4, 'filename': 'orbital_4.cube'},
                {'orbital_index': 5, 'filename': 'orbital_5.cube'}
            ],
            'total_files': 2
        }
        mock_service = mocker.patch('api.quantum.get_quantum_service')
        mock_service.return_value.list_cube_files.return_value = mock_files

        # ACT
        response = client.get(f'/api/quantum/calculations/{calc_id}/orbitals/cube-files')

        # ASSERT
        assert response.status_code == 200
        data = response.get_json()
        assert data['success'] is True
        assert data['data']['total_files'] == 2

    def test_delete_cube_files_all(self, client, mocker):
        """
        GIVEN calculation has CUBE files
        WHEN DELETE /api/quantum/calculations/<id>/orbitals/cube-files is called
        THEN all CUBE files are deleted
        """
        # ARRANGE
        calc_id = 'calc-123'
        mock_result = {
            'deleted_files': 5,
            'message': 'All CUBE files deleted'
        }
        mock_service = mocker.patch('api.quantum.get_quantum_service')
        mock_service.return_value.delete_cube_files.return_value = mock_result

        # ACT
        response = client.delete(f'/api/quantum/calculations/{calc_id}/orbitals/cube-files')

        # ASSERT
        assert response.status_code == 200
        data = response.get_json()
        assert data['success'] is True
        assert data['data']['deleted_files'] == 5

    def test_delete_cube_files_specific_orbital(self, client, mocker):
        """
        GIVEN orbital_index query parameter
        WHEN DELETE is called
        THEN only that orbital's CUBE file is deleted
        """
        # ARRANGE
        calc_id = 'calc-123'
        orbital_index = 5
        mock_result = {
            'deleted_files': 1,
            'orbital_index': orbital_index
        }
        mock_service = mocker.patch('api.quantum.get_quantum_service')
        mock_service.return_value.delete_cube_files.return_value = mock_result

        # ACT
        response = client.delete(
            f'/api/quantum/calculations/{calc_id}/orbitals/cube-files'
            f'?orbital_index={orbital_index}'
        )

        # ASSERT
        assert response.status_code == 200
        mock_service.return_value.delete_cube_files.assert_called_once_with(calc_id, orbital_index)


class TestIRSpectrumAPI:
    """Integration tests for IR spectrum generation endpoint."""

    def test_generate_ir_spectrum_success(self, client, mocker):
        """
        GIVEN calculation has frequency data
        WHEN GET /api/quantum/calculations/<id>/ir-spectrum is called
        THEN 200 OK is returned with spectrum data
        """
        # ARRANGE
        calc_id = 'calc-123'
        mock_spectrum = {
            'spectrum': {
                'x': [400, 500, 600],
                'y': [0.1, 0.5, 0.2]
            },
            'plot_image_base64': 'base64encodedimage...'
        }
        mock_service = mocker.patch('api.quantum.get_quantum_service')
        mock_service.return_value.generate_ir_spectrum.return_value = mock_spectrum

        # ACT
        response = client.get(f'/api/quantum/calculations/{calc_id}/ir-spectrum')

        # ASSERT
        assert response.status_code == 200
        data = response.get_json()
        assert data['success'] is True
        assert 'spectrum' in data['data']

    def test_generate_ir_spectrum_with_custom_params(self, client, mocker):
        """
        GIVEN custom spectrum generation parameters
        WHEN GET with query parameters is called
        THEN service receives the custom parameters
        """
        # ARRANGE
        calc_id = 'calc-123'
        mock_spectrum = {'spectrum': {}}
        mock_service = mocker.patch('api.quantum.get_quantum_service')
        mock_service.return_value.generate_ir_spectrum.return_value = mock_spectrum

        # ACT
        response = client.get(
            f'/api/quantum/calculations/{calc_id}/ir-spectrum'
            f'?broadening_fwhm=50&x_min=500&x_max=3500'
        )

        # ASSERT
        assert response.status_code == 200
        mock_service.return_value.generate_ir_spectrum.assert_called_once_with(
            calc_id,
            broadening_fwhm=50.0,
            x_min=500.0,
            x_max=3500.0,
            show_peaks=True  # Default value when not specified
        )

    def test_generate_ir_spectrum_not_found(self, client, mocker):
        """
        GIVEN calculation does not have frequency data
        WHEN GET /api/quantum/calculations/<id>/ir-spectrum is called
        THEN 404 Not Found is returned
        """
        # ARRANGE
        calc_id = 'calc-no-freq'
        mock_service = mocker.patch('api.quantum.get_quantum_service')
        mock_service.return_value.generate_ir_spectrum.side_effect = NotFoundError(f"No frequency data available")

        # ACT
        response = client.get(f'/api/quantum/calculations/{calc_id}/ir-spectrum')

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
            'process_pool': {
                'active': True,
                'workers': 4
            },
            'system': {
                'cpu_count': 8,
                'memory_available': 16000
            }
        }
        mock_service = mocker.patch('api.quantum.get_quantum_service')
        mock_service.return_value.get_calculation_status.return_value = mock_status

        # ACT
        response = client.get('/api/quantum/status')

        # ASSERT
        assert response.status_code == 200
        data = response.get_json()
        assert data['success'] is True
        assert 'process_pool' in data['data']
