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
            'xyz': 'O 0.0000 0.0000 0.1173\nH 0.0000 0.7572 -0.4692\nH 0.0000 -0.7572 -0.4692',
            'cid': 962,
            'name': 'Water'
        }
        mock_service = mocker.patch('api.pubchem.get_pubchem_service')
        mock_service.return_value.search_compound.return_value = mock_result

        # ACT
        response = client.post('/api/pubchem/search', json={
            'query': 'water',
            'searchType': 'name'
        })

        # ASSERT
        assert response.status_code == 200
        data = response.get_json()
        assert data['success'] is True
        assert 'data' in data
        assert data['data']['xyz'] == mock_result['xyz']
        assert data['data']['cid'] == 962
        
        # Verify service was called correctly
        mock_service.return_value.search_compound.assert_called_once_with('water', 'name')

    @pytest.mark.parametrize("search_type", ['name', 'cid', 'formula'])
    def test_search_different_types(self, client, mocker, search_type):
        """
        GIVEN PubChemService is configured
        WHEN POST /api/pubchem/search is called with different search types
        THEN the service is called with the correct search type
        """
        # ARRANGE
        mock_result = {'xyz': 'H 0 0 0', 'cid': 123}
        mock_service = mocker.patch('api.pubchem.get_pubchem_service')
        mock_service.return_value.search_compound.return_value = mock_result

        # ACT
        response = client.post('/api/pubchem/search', json={
            'query': 'test_query',
            'searchType': search_type
        })

        # ASSERT
        assert response.status_code == 200
        mock_service.return_value.search_compound.assert_called_once_with('test_query', search_type)

    def test_search_not_found(self, client, mocker):
        """
        GIVEN PubChemService raises NotFoundError
        WHEN POST /api/pubchem/search is called
        THEN 404 Not Found is returned
        """
        # ARRANGE
        mock_service = mocker.patch('api.pubchem.get_pubchem_service')
        mock_service.return_value.search_compound.side_effect = NotFoundError("Compound not found")

        # ACT
        response = client.post('/api/pubchem/search', json={
            'query': 'nonexistent_compound_xyz123',
            'searchType': 'name'
        })

        # ASSERT
        assert response.status_code == 404
        data = response.get_json()
        assert data['success'] is False
        assert 'error' in data
        assert 'not found' in data['error'].lower()

    def test_search_service_error(self, client, mocker):
        """
        GIVEN PubChemService raises ServiceError
        WHEN POST /api/pubchem/search is called
        THEN appropriate error status is returned
        """
        # ARRANGE
        mock_service = mocker.patch('api.pubchem.get_pubchem_service')
        mock_service.return_value.search_compound.side_effect = ServiceError(
            "PubChem API unavailable",
            status_code=503
        )

        # ACT
        response = client.post('/api/pubchem/search', json={
            'query': 'water',
            'searchType': 'name'
        })

        # ASSERT
        assert response.status_code == 503
        data = response.get_json()
        assert data['success'] is False
        assert 'error' in data

    def test_search_missing_required_fields(self, client):
        """
        GIVEN request is missing required fields
        WHEN POST /api/pubchem/search is called
        THEN 400 Bad Request is returned
        """
        # ACT - Missing query (required field)
        response = client.post('/api/pubchem/search', json={
            'searchType': 'name'
        })

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
            '/api/pubchem/search',
            data='invalid json',
            content_type='application/json'
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
        valid_xyz = 'H 0 0 0\nH 0 0 0.74'
        mock_result = {
            'valid': True,
            'atom_count': 2,
            'message': 'XYZ format is valid'
        }
        mock_service = mocker.patch('api.pubchem.get_pubchem_service')
        mock_service.return_value.validate_xyz.return_value = mock_result

        # ACT
        response = client.post('/api/pubchem/validate', json={
            'xyz': valid_xyz
        })

        # ASSERT
        assert response.status_code == 200
        data = response.get_json()
        assert data['success'] is True
        assert data['data']['valid'] is True
        assert data['data']['atom_count'] == 2
        
        mock_service.return_value.validate_xyz.assert_called_once_with(valid_xyz)

    def test_validate_invalid_xyz(self, client, mocker):
        """
        GIVEN PubChemService returns validation failure
        WHEN POST /api/pubchem/validate is called with invalid XYZ
        THEN 200 OK is returned with validation error details
        """
        # ARRANGE
        invalid_xyz = 'invalid xyz format'
        mock_result = {
            'valid': False,
            'message': 'Invalid XYZ format'
        }
        mock_service = mocker.patch('api.pubchem.get_pubchem_service')
        mock_service.return_value.validate_xyz.return_value = mock_result

        # ACT
        response = client.post('/api/pubchem/validate', json={
            'xyz': invalid_xyz
        })

        # ASSERT
        assert response.status_code == 200
        data = response.get_json()
        assert data['success'] is True
        assert data['data']['valid'] is False

    def test_validate_empty_xyz(self, client):
        """
        GIVEN empty XYZ string
        WHEN POST /api/pubchem/validate is called
        THEN 400 Bad Request is returned
        """
        # ACT
        response = client.post('/api/pubchem/validate', json={
            'xyz': ''
        })

        # ASSERT
        assert response.status_code == 400

    def test_validate_missing_xyz_field(self, client):
        """
        GIVEN request is missing xyz field
        WHEN POST /api/pubchem/validate is called
        THEN 400 Bad Request is returned
        """
        # ACT
        response = client.post('/api/pubchem/validate', json={})

        # ASSERT
        assert response.status_code == 400
