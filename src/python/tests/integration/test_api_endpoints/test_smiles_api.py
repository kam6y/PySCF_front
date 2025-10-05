"""
Integration tests for SMILES Conversion API endpoints.

Tests the SMILES-to-XYZ conversion endpoint, ensuring proper
handling of various SMILES strings and error conditions.
"""

import pytest
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
            'xyz': 'O 0.0000 0.0000 0.1173\nH 0.0000 0.7572 -0.4692\nH 0.0000 -0.7572 -0.4692',
            'smiles': 'O'
        }
        mock_service = mocker.patch('api.smiles.get_smiles_service')
        mock_service.return_value.convert_smiles.return_value = mock_result

        # ACT
        response = client.post('/api/smiles/convert', json={
            'smiles': 'O'
        })

        # ASSERT
        assert response.status_code == 200
        data = response.get_json()
        assert data['success'] is True
        assert 'data' in data
        assert 'xyz' in data['data']
        assert 'O' in data['data']['xyz']
        assert 'H' in data['data']['xyz']
        
        # Verify service was called correctly
        mock_service.return_value.convert_smiles.assert_called_once_with('O')

    def test_convert_complex_smiles(self, client, mocker):
        """
        GIVEN SMILESService handles complex molecules
        WHEN POST /api/smiles/convert is called with benzene SMILES
        THEN conversion succeeds
        """
        # ARRANGE
        benzene_smiles = 'c1ccccc1'
        mock_result = {
            'xyz': 'C 0 0 0\nC 1 0 0\nC 1.5 0.866 0\nC 1.5 1.732 0\nC 1 2.598 0\nC 0 2.598 0',
            'smiles': benzene_smiles
        }
        mock_service = mocker.patch('api.smiles.get_smiles_service')
        mock_service.return_value.convert_smiles.return_value = mock_result

        # ACT
        response = client.post('/api/smiles/convert', json={
            'smiles': benzene_smiles
        })

        # ASSERT
        assert response.status_code == 200
        data = response.get_json()
        assert data['success'] is True
        assert 'C' in data['data']['xyz']

    def test_convert_with_whitespace_trimming(self, client, mocker):
        """
        GIVEN SMILES string has leading/trailing whitespace
        WHEN POST /api/smiles/convert is called
        THEN whitespace is trimmed before conversion
        """
        # ARRANGE
        smiles_with_whitespace = '  CCO  '
        mock_result = {'xyz': 'C 0 0 0\nC 1 0 0\nO 2 0 0', 'smiles': 'CCO'}
        mock_service = mocker.patch('api.smiles.get_smiles_service')
        mock_service.return_value.convert_smiles.return_value = mock_result

        # ACT
        response = client.post('/api/smiles/convert', json={
            'smiles': smiles_with_whitespace
        })

        # ASSERT
        assert response.status_code == 200
        # Service should receive trimmed version
        mock_service.return_value.convert_smiles.assert_called_once_with(smiles_with_whitespace)

    def test_convert_invalid_smiles(self, client, mocker):
        """
        GIVEN SMILESService raises ValidationError for invalid SMILES
        WHEN POST /api/smiles/convert is called
        THEN 400 Bad Request is returned
        """
        # ARRANGE
        invalid_smiles = 'INVALID_SMILES_XXX'
        mock_service = mocker.patch('api.smiles.get_smiles_service')
        mock_service.return_value.convert_smiles.side_effect = ValidationError("Invalid SMILES string")

        # ACT
        response = client.post('/api/smiles/convert', json={
            'smiles': invalid_smiles
        })

        # ASSERT
        assert response.status_code == 400
        data = response.get_json()
        assert data['success'] is False
        assert 'error' in data

    def test_convert_empty_smiles(self, client):
        """
        GIVEN empty SMILES string
        WHEN POST /api/smiles/convert is called
        THEN 400 Bad Request is returned
        """
        # ACT
        response = client.post('/api/smiles/convert', json={
            'smiles': ''
        })

        # ASSERT
        assert response.status_code == 400

    def test_convert_missing_smiles_field(self, client):
        """
        GIVEN request is missing smiles field
        WHEN POST /api/smiles/convert is called
        THEN 400 Bad Request is returned
        """
        # ACT
        response = client.post('/api/smiles/convert', json={})

        # ASSERT
        assert response.status_code == 400

    def test_convert_service_error(self, client, mocker):
        """
        GIVEN SMILESService raises ServiceError
        WHEN POST /api/smiles/convert is called
        THEN appropriate error status is returned
        """
        # ARRANGE
        mock_service = mocker.patch('api.smiles.get_smiles_service')
        mock_service.return_value.convert_smiles.side_effect = ServiceError(
            "RDKit library error",
            status_code=500
        )

        # ACT
        response = client.post('/api/smiles/convert', json={
            'smiles': 'CCO'
        })

        # ASSERT
        assert response.status_code == 500
        data = response.get_json()
        assert data['success'] is False

    def test_convert_whitespace_only_smiles(self, client):
        """
        GIVEN SMILES string contains only whitespace
        WHEN POST /api/smiles/convert is called
        THEN 400 Bad Request is returned
        """
        # ACT
        response = client.post('/api/smiles/convert', json={
            'smiles': '   '
        })

        # ASSERT
        assert response.status_code == 400

    def test_convert_very_long_smiles(self, client, mocker):
        """
        GIVEN very long SMILES string
        WHEN POST /api/smiles/convert is called
        THEN service handles it appropriately
        """
        # ARRANGE
        # Create a long SMILES (polymer-like)
        long_smiles = 'C' * 1000
        mock_result = {'xyz': 'C 0 0 0', 'smiles': long_smiles}
        mock_service = mocker.patch('api.smiles.get_smiles_service')
        mock_service.return_value.convert_smiles.return_value = mock_result

        # ACT
        response = client.post('/api/smiles/convert', json={
            'smiles': long_smiles
        })

        # ASSERT
        assert response.status_code == 200
        data = response.get_json()
        assert data['success'] is True

    def test_convert_invalid_json(self, client):
        """
        GIVEN invalid JSON payload
        WHEN POST /api/smiles/convert is called
        THEN 400 Bad Request is returned
        """
        # ACT
        response = client.post(
            '/api/smiles/convert',
            data='invalid json',
            content_type='application/json'
        )

        # ASSERT
        assert response.status_code == 400
