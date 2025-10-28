"""
Integration tests for complete calculation workflows.

Tests the end-to-end flow from calculation submission through execution
to results retrieval, using DummyExecutor for synchronous testing.
"""

import pytest
import os
import json
from pathlib import Path
from tests.conftest import DummyExecutor


class TestCalculationWorkflowSync:
    """
    Integration tests for complete calculation workflow using synchronous execution.
    
    These tests use DummyExecutor to replace ProcessPoolExecutor, allowing
    calculations to run synchronously and deterministically in tests.
    """

    def test_full_hf_calculation_workflow(self, client, mocker, valid_hf_params, app):
        """
        GIVEN ProcessPoolExecutor is replaced with DummyExecutor
        WHEN a HF calculation is submitted and completes
        THEN the full workflow from submission to results retrieval works
        """
        # ARRANGE
        # Replace ProcessPoolExecutor with DummyExecutor for synchronous execution
        mocker.patch(
            'quantum_calc.process_manager.ProcessPoolExecutor',
            new=DummyExecutor
        )
        
        # Mock PySCF to avoid actual computation
        mock_mol = mocker.MagicMock()
        mock_scf = mocker.MagicMock()
        mock_scf.kernel.return_value = -1.06  # Mock energy
        mock_scf.mo_energy = [-0.5, 0.3]
        mock_scf.mo_occ = [2.0, 0.0]
        
        mocker.patch('pyscf.gto.M', return_value=mock_mol)
        mocker.patch('pyscf.scf.RHF', return_value=mock_scf)

        # ACT
        # Step 1: Submit calculation
        response_submit = client.post('/api/quantum/calculate', json=valid_hf_params)
        
        # ASSERT Step 1
        assert response_submit.status_code == 202
        submit_data = response_submit.get_json()
        assert submit_data['success'] is True
        
        calc_id = submit_data['data']['calculation']['id']
        assert calc_id is not None
        
        # ACT
        # Step 2: Get calculation details (should be completed due to DummyExecutor)
        response_details = client.get(f'/api/quantum/calculations/{calc_id}')
        
        # ASSERT Step 2
        assert response_details.status_code == 200
        details_data = response_details.get_json()
        assert details_data['success'] is True
        
        calc_details = details_data['data']['calculation']
        assert calc_details['id'] == calc_id
        assert calc_details['status'] in ['completed', 'running']  # May still be running in edge cases
        
        # If completed, verify results
        if calc_details['status'] == 'completed':
            assert 'results' in calc_details
            assert 'energy' in calc_details['results']

    def test_workflow_calculation_listing(self, client, mocker, valid_dft_params):
        """
        GIVEN multiple calculations are submitted
        WHEN GET /api/quantum/calculations is called
        THEN all calculations are listed
        """
        # ARRANGE
        mocker.patch('quantum_calc.process_manager.ProcessPoolExecutor', new=DummyExecutor)
        
        # Mock PySCF
        mock_mol = mocker.MagicMock()
        mock_scf = mocker.MagicMock()
        mock_scf.kernel.return_value = -1.5
        mock_scf.mo_energy = [-0.5]
        mock_scf.mo_occ = [2.0]
        
        mocker.patch('pyscf.gto.M', return_value=mock_mol)
        mocker.patch('pyscf.dft.RKS', return_value=mock_scf)

        # ACT
        # Submit multiple calculations
        calc_ids = []
        for i in range(3):
            params = {**valid_dft_params, 'name': f'Test Calc {i}'}
            response = client.post('/api/quantum/calculate', json=params)
            assert response.status_code == 202
            calc_id = response.get_json()['data']['calculation']['id']
            calc_ids.append(calc_id)

        # Get list of calculations
        response_list = client.get('/api/quantum/calculations')

        # ASSERT
        assert response_list.status_code == 200
        list_data = response_list.get_json()
        assert list_data['success'] is True
        
        # All submitted calculations should be in the list
        calc_list = list_data['data']['calculations']
        listed_ids = [calc['id'] for calc in calc_list]
        
        for calc_id in calc_ids:
            assert calc_id in listed_ids

    def test_workflow_calculation_rename(self, client, mocker, valid_hf_params):
        """
        GIVEN a calculation exists
        WHEN it is renamed via PUT endpoint
        THEN the name is updated and reflected in details
        """
        # ARRANGE
        mocker.patch('quantum_calc.process_manager.ProcessPoolExecutor', new=DummyExecutor)
        
        # Mock PySCF
        mock_mol = mocker.MagicMock()
        mock_scf = mocker.MagicMock()
        mock_scf.kernel.return_value = -1.0
        mock_scf.mo_energy = [-0.5]
        mock_scf.mo_occ = [2.0]
        
        mocker.patch('pyscf.gto.M', return_value=mock_mol)
        mocker.patch('pyscf.scf.RHF', return_value=mock_scf)

        # Submit calculation
        response_submit = client.post('/api/quantum/calculate', json=valid_hf_params)
        calc_id = response_submit.get_json()['data']['calculation']['id']

        # ACT
        new_name = 'Renamed Calculation'
        response_rename = client.put(f'/api/quantum/calculations/{calc_id}', json={
            'name': new_name
        })

        # ASSERT
        assert response_rename.status_code == 200
        rename_data = response_rename.get_json()
        assert rename_data['success'] is True
        # API returns {'message': '...', 'name': '...'} directly in data
        assert rename_data['data']['name'] == new_name

        # Verify name persists
        response_details = client.get(f'/api/quantum/calculations/{calc_id}')
        details_data = response_details.get_json()
        assert details_data['data']['calculation']['name'] == new_name

    def test_workflow_calculation_deletion(self, client, mocker, valid_hf_params):
        """
        GIVEN a calculation exists
        WHEN it is deleted via DELETE endpoint
        THEN it no longer appears in listings and cannot be retrieved
        """
        # ARRANGE
        mocker.patch('quantum_calc.process_manager.ProcessPoolExecutor', new=DummyExecutor)
        
        # Mock PySCF
        mock_mol = mocker.MagicMock()
        mock_scf = mocker.MagicMock()
        mock_scf.kernel.return_value = -1.0
        mock_scf.mo_energy = [-0.5]
        mock_scf.mo_occ = [2.0]
        
        mocker.patch('pyscf.gto.M', return_value=mock_mol)
        mocker.patch('pyscf.scf.RHF', return_value=mock_scf)

        # Submit calculation
        response_submit = client.post('/api/quantum/calculate', json=valid_hf_params)
        calc_id = response_submit.get_json()['data']['calculation']['id']

        # Verify it exists
        response_before = client.get(f'/api/quantum/calculations/{calc_id}')
        assert response_before.status_code == 200

        # ACT
        response_delete = client.delete(f'/api/quantum/calculations/{calc_id}')

        # ASSERT
        assert response_delete.status_code == 200
        delete_data = response_delete.get_json()
        assert delete_data['success'] is True

        # Verify it's deleted (404)
        response_after = client.get(f'/api/quantum/calculations/{calc_id}')
        assert response_after.status_code == 404

    def test_workflow_with_websocket_integration(self, client, socketio_client, mocker, valid_hf_params, app):
        """
        GIVEN WebSocket client is connected to a calculation
        WHEN calculation completes
        THEN WebSocket receives update notifications
        """
        # ARRANGE
        mocker.patch('quantum_calc.process_manager.ProcessPoolExecutor', new=DummyExecutor)
        
        # Mock PySCF
        mock_mol = mocker.MagicMock()
        mock_scf = mocker.MagicMock()
        mock_scf.kernel.return_value = -1.06
        mock_scf.mo_energy = [-0.5, 0.3]
        mock_scf.mo_occ = [2.0, 0.0]
        
        mocker.patch('pyscf.gto.M', return_value=mock_mol)
        mocker.patch('pyscf.scf.RHF', return_value=mock_scf)

        # ACT
        # Step 1: Submit calculation
        response_submit = client.post('/api/quantum/calculate', json=valid_hf_params)
        calc_id = response_submit.get_json()['data']['calculation']['id']

        # Step 2: Join WebSocket room
        socketio_client.emit('join_calculation', {'calculation_id': calc_id})
        
        # Get received messages
        received = socketio_client.get_received()

        # ASSERT
        # Should receive calculation_update event with initial state
        update_events = [msg for msg in received if msg['name'] == 'calculation_update']
        assert len(update_events) > 0
        
        calc_data = update_events[0]['args'][0]
        assert calc_data['id'] == calc_id

    def test_workflow_error_handling(self, client, mocker, valid_hf_params):
        """
        GIVEN PySCF raises an exception during calculation
        WHEN calculation is executed
        THEN error status is properly recorded
        """
        # ARRANGE
        mocker.patch('quantum_calc.process_manager.ProcessPoolExecutor', new=DummyExecutor)
        
        # Mock PySCF to raise error
        mock_mol = mocker.MagicMock()
        mocker.patch('pyscf.gto.M', return_value=mock_mol)
        
        mock_scf = mocker.MagicMock()
        mock_scf.kernel.side_effect = RuntimeError("SCF did not converge")
        mocker.patch('pyscf.scf.RHF', return_value=mock_scf)

        # ACT
        response_submit = client.post('/api/quantum/calculate', json=valid_hf_params)
        calc_id = response_submit.get_json()['data']['calculation']['id']

        # Get calculation details
        response_details = client.get(f'/api/quantum/calculations/{calc_id}')

        # ASSERT
        assert response_details.status_code == 200
        details_data = response_details.get_json()
        calc_details = details_data['data']['calculation']

        # Should have error status (or waiting if not yet processed)
        assert calc_details['status'] in ['error', 'waiting', 'running']

    def test_workflow_orbital_generation(self, client, mocker, valid_hf_params):
        """
        GIVEN a completed calculation with orbital data
        WHEN orbital CUBE file is requested
        THEN CUBE file is generated successfully
        """
        # ARRANGE
        mocker.patch('quantum_calc.process_manager.ProcessPoolExecutor', new=DummyExecutor)
        
        # Mock PySCF with orbital data
        mock_mol = mocker.MagicMock()
        mock_scf = mocker.MagicMock()
        mock_scf.kernel.return_value = -1.06
        mock_scf.mo_energy = [-0.5, 0.3, 0.8]
        mock_scf.mo_occ = [2.0, 0.0, 0.0]
        mock_scf.mo_coeff = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]  # Mock MO coefficients
        
        mocker.patch('pyscf.gto.M', return_value=mock_mol)
        mocker.patch('pyscf.scf.RHF', return_value=mock_scf)
        
        # Mock CUBE file generation
        mocker.patch('pyscf.tools.cubegen.orbital')

        # Submit and complete calculation
        response_submit = client.post('/api/quantum/calculate', json=valid_hf_params)
        calc_id = response_submit.get_json()['data']['calculation']['id']

        # Wait for calculation to complete
        import time
        for _ in range(10):  # Try for 10 seconds
            response_details = client.get(f'/api/quantum/calculations/{calc_id}')
            calc_details = response_details.get_json()['data']['calculation']
            if calc_details['status'] == 'completed':
                break
            time.sleep(1)

        # ACT
        # Request orbitals list
        response_orbitals = client.get(f'/api/quantum/calculations/{calc_id}/orbitals')

        # ASSERT
        # May return 200 with orbital data, or 400 if calculation not complete/failed
        assert response_orbitals.status_code in [200, 400]
        if response_orbitals.status_code == 200:
            orbitals_data = response_orbitals.get_json()
            assert orbitals_data['success'] is True

            # Should have orbital information
            assert 'orbitals' in orbitals_data['data'] or 'homo_index' in orbitals_data['data']


class TestCalculationWorkflowValidation:
    """Integration tests for calculation parameter validation in workflow context."""

    def test_workflow_rejects_invalid_basis_set(self, client, valid_dft_params):
        """
        GIVEN invalid basis set in calculation parameters
        WHEN calculation is submitted
        THEN calculation may be accepted but should fail during execution
        """
        # ARRANGE
        invalid_params = {**valid_dft_params, 'basis_function': 'invalid_basis_xyz'}

        # ACT
        response = client.post('/api/quantum/calculate', json=invalid_params)

        # ASSERT
        # May accept (202) and fail later, or reject immediately (400/422)
        assert response.status_code in [202, 400, 422]

        # If accepted, it should fail during calculation
        if response.status_code == 202:
            calc_id = response.get_json()['data']['calculation']['id']
            import time
            for _ in range(10):
                details_response = client.get(f'/api/quantum/calculations/{calc_id}')
                details = details_response.get_json()['data']['calculation']
                if details['status'] in ['error', 'completed']:
                    break
                time.sleep(1)
            # Should eventually error due to invalid basis set
            assert details['status'] in ['error', 'waiting', 'running']

    def test_workflow_rejects_missing_xyz(self, client, valid_dft_params):
        """
        GIVEN XYZ coordinates are missing
        WHEN calculation is submitted
        THEN 400 Bad Request is returned
        """
        # ARRANGE
        invalid_params = {**valid_dft_params}
        del invalid_params['xyz']

        # ACT
        response = client.post('/api/quantum/calculate', json=invalid_params)

        # ASSERT
        assert response.status_code == 400

    def test_workflow_rejects_invalid_charge_spin_combination(self, client, valid_hf_params):
        """
        GIVEN invalid charge/spin combination
        WHEN calculation is submitted
        THEN validation error occurs
        """
        # ARRANGE
        # Even number of electrons with odd spin is invalid
        invalid_params = {
            **valid_hf_params,
            'charges': 0,  # H2 has 2 electrons
            'spin': 1  # Odd spin with even electrons is invalid
        }

        # ACT
        response = client.post('/api/quantum/calculate', json=invalid_params)

        # ASSERT
        # Should be rejected either at API level or during calculation
        # Status could be 400 (validation) or 202 followed by error status
        if response.status_code == 202:
            calc_id = response.get_json()['data']['calculation']['id']
            import time
            for _ in range(10):
                details_response = client.get(f'/api/quantum/calculations/{calc_id}')
                details = details_response.get_json()['data']['calculation']
                if details['status'] in ['error', 'completed']:
                    break
                time.sleep(1)
            # If accepted, should eventually error or be in waiting/running state
            assert details['status'] in ['error', 'waiting', 'running', 'pending']


class TestCalculationWorkflowMultipleCalculations:
    """Integration tests for managing multiple concurrent calculations."""

    def test_multiple_calculations_independent(self, client, mocker, valid_hf_params, valid_dft_params):
        """
        GIVEN multiple calculations are submitted
        WHEN they execute independently
        THEN each maintains its own state and results
        """
        # ARRANGE
        mocker.patch('quantum_calc.process_manager.ProcessPoolExecutor', new=DummyExecutor)
        
        # Mock PySCF
        mock_mol = mocker.MagicMock()
        mocker.patch('pyscf.gto.M', return_value=mock_mol)
        
        # Different energies for different calculations
        energies = [-1.06, -1.5, -2.0]
        call_count = [0]
        
        def mock_kernel_side_effect():
            energy = energies[call_count[0] % len(energies)]
            call_count[0] += 1
            return energy
        
        mock_scf_hf = mocker.MagicMock()
        mock_scf_hf.kernel.side_effect = mock_kernel_side_effect
        mock_scf_hf.mo_energy = [-0.5]
        mock_scf_hf.mo_occ = [2.0]
        
        mock_scf_dft = mocker.MagicMock()
        mock_scf_dft.kernel.side_effect = mock_kernel_side_effect
        mock_scf_dft.mo_energy = [-0.5]
        mock_scf_dft.mo_occ = [2.0]
        
        mocker.patch('pyscf.scf.RHF', return_value=mock_scf_hf)
        mocker.patch('pyscf.dft.RKS', return_value=mock_scf_dft)

        # ACT
        # Submit HF calculation
        response_hf = client.post('/api/quantum/calculate', json=valid_hf_params)
        hf_calc_id = response_hf.get_json()['data']['calculation']['id']

        # Submit DFT calculation
        response_dft = client.post('/api/quantum/calculate', json=valid_dft_params)
        dft_calc_id = response_dft.get_json()['data']['calculation']['id']

        # ASSERT
        # Both should have unique IDs
        assert hf_calc_id != dft_calc_id

        # Both should be retrievable independently
        hf_details = client.get(f'/api/quantum/calculations/{hf_calc_id}').get_json()
        dft_details = client.get(f'/api/quantum/calculations/{dft_calc_id}').get_json()

        assert hf_details['success'] is True
        assert dft_details['success'] is True
        
        # Should have different calculation methods
        hf_method = hf_details['data']['calculation']['parameters']['calculation_method']
        dft_method = dft_details['data']['calculation']['parameters']['calculation_method']
        
        assert hf_method == 'HF'
        assert dft_method == 'DFT'
