"""Tests for Flask API endpoints."""

import pytest
import json
from unittest.mock import Mock, patch

# Add the parent directory to the path for imports
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from app import app
from pubchem.client import PubChemError


class TestFlaskAPI:
    """Test cases for Flask API endpoints."""
    
    @pytest.fixture
    def client(self):
        """Create test client."""
        app.config['TESTING'] = True
        with app.test_client() as client:
            yield client
    
    def test_health_check(self, client):
        """Test health check endpoint."""
        response = client.get('/health')
        assert response.status_code == 200
        
        data = json.loads(response.data)
        assert data['status'] == 'ok'
        assert data['service'] == 'pyscf-pubchem-api'
    
    @patch('app.pubchem_client.search_compound')
    @patch('app.pubchem_client.get_compound_info')
    @patch('app.pubchem_client.get_3d_structure')
    @patch('app.xyz_parser.atoms_to_xyz')
    @patch('app.xyz_parser.validate_xyz')
    def test_search_pubchem_success(self, mock_validate, mock_to_xyz, 
                                  mock_get_3d, mock_get_info, mock_search, client):
        """Test successful PubChem search."""
        # Mock compound
        mock_compound = Mock()
        mock_compound.cid = 962
        mock_search.return_value = mock_compound
        
        # Mock compound info
        mock_info = {
            'cid': 962,
            'molecular_formula': 'H2O',
            'iupac_name': 'water'
        }
        mock_get_info.return_value = mock_info
        
        # Mock 3D structure
        mock_atoms = [
            {'element': 'O', 'x': 0.0, 'y': 0.0, 'z': 0.1193}
        ]
        mock_get_3d.return_value = mock_atoms
        
        # Mock XYZ conversion
        mock_xyz = "1\nwater\nO 0.000000 0.000000 0.119300"
        mock_to_xyz.return_value = mock_xyz
        
        # Mock XYZ validation
        mock_validate.return_value = {'valid': True}
        
        # Make request
        response = client.post('/api/pubchem/search',
                             json={'query': 'water', 'search_type': 'name'},
                             content_type='application/json')
        
        assert response.status_code == 200
        
        data = json.loads(response.data)
        assert data['success'] is True
        assert data['data']['xyz'] == mock_xyz
        assert data['data']['compound_info'] == mock_info
        assert data['data']['atom_count'] == 1
    
    def test_search_pubchem_no_json(self, client):
        """Test PubChem search with no JSON data."""
        response = client.post('/api/pubchem/search')
        
        assert response.status_code == 400
        
        data = json.loads(response.data)
        assert data['success'] is False
        assert 'No JSON data provided' in data['error']
    
    def test_search_pubchem_empty_query(self, client):
        """Test PubChem search with empty query."""
        response = client.post('/api/pubchem/search',
                             json={'query': ''},
                             content_type='application/json')
        
        assert response.status_code == 400
        
        data = json.loads(response.data)
        assert data['success'] is False
        assert 'Query parameter is required' in data['error']
    
    def test_search_pubchem_invalid_search_type(self, client):
        """Test PubChem search with invalid search type."""
        response = client.post('/api/pubchem/search',
                             json={'query': 'water', 'search_type': 'invalid'},
                             content_type='application/json')
        
        assert response.status_code == 400
        
        data = json.loads(response.data)
        assert data['success'] is False
        assert 'Invalid search_type' in data['error']
    
    @patch('app.pubchem_client.search_compound')
    def test_search_pubchem_compound_not_found(self, mock_search, client):
        """Test PubChem search when compound not found."""
        mock_search.return_value = None
        
        response = client.post('/api/pubchem/search',
                             json={'query': 'nonexistent'},
                             content_type='application/json')
        
        assert response.status_code == 404
        
        data = json.loads(response.data)
        assert data['success'] is False
        assert 'No compound found' in data['error']
    
    @patch('app.pubchem_client.search_compound')
    @patch('app.pubchem_client.get_3d_structure')
    def test_search_pubchem_no_3d_structure(self, mock_get_3d, mock_search, client):
        """Test PubChem search when no 3D structure available."""
        mock_compound = Mock()
        mock_compound.cid = 123
        mock_search.return_value = mock_compound
        mock_get_3d.return_value = None
        
        response = client.post('/api/pubchem/search',
                             json={'query': 'test'},
                             content_type='application/json')
        
        assert response.status_code == 404
        
        data = json.loads(response.data)
        assert data['success'] is False
        assert 'No 3D structure available' in data['error']
    
    @patch('app.pubchem_client.search_compound')
    def test_search_pubchem_api_error(self, mock_search, client):
        """Test PubChem search with API error."""
        mock_search.side_effect = PubChemError("API error")
        
        response = client.post('/api/pubchem/search',
                             json={'query': 'water'},
                             content_type='application/json')
        
        assert response.status_code == 500
        
        data = json.loads(response.data)
        assert data['success'] is False
        assert 'API error' in data['error']
    
    @patch('app.xyz_parser.validate_xyz')
    def test_validate_xyz_endpoint_success(self, mock_validate, client):
        """Test XYZ validation endpoint success."""
        mock_validate.return_value = {
            'valid': True,
            'num_atoms': 3,
            'title': 'Water'
        }
        
        xyz_string = "3\nWater\nO 0.0 0.0 0.0"
        response = client.post('/api/pubchem/validate',
                             json={'xyz': xyz_string},
                             content_type='application/json')
        
        assert response.status_code == 200
        
        data = json.loads(response.data)
        assert data['success'] is True
        assert data['data']['valid'] is True
        assert data['data']['num_atoms'] == 3
    
    @patch('app.xyz_parser.validate_xyz')
    def test_validate_xyz_endpoint_invalid(self, mock_validate, client):
        """Test XYZ validation endpoint with invalid XYZ."""
        mock_validate.return_value = {
            'valid': False,
            'error': 'Invalid format'
        }
        
        response = client.post('/api/pubchem/validate',
                             json={'xyz': 'invalid'},
                             content_type='application/json')
        
        assert response.status_code == 200
        
        data = json.loads(response.data)
        assert data['success'] is True
        assert data['data']['valid'] is False
        assert 'Invalid format' in data['data']['error']
    
    def test_validate_xyz_no_data(self, client):
        """Test XYZ validation with no JSON data."""
        response = client.post('/api/pubchem/validate')
        
        assert response.status_code == 400
        
        data = json.loads(response.data)
        assert data['success'] is False
        assert 'No JSON data provided' in data['error']
    
    def test_validate_xyz_empty_string(self, client):
        """Test XYZ validation with empty string."""
        response = client.post('/api/pubchem/validate',
                             json={'xyz': ''},
                             content_type='application/json')
        
        assert response.status_code == 400
        
        data = json.loads(response.data)
        assert data['success'] is False
        assert 'XYZ string is required' in data['error']
    
    def test_not_found_endpoint(self, client):
        """Test 404 error handler."""
        response = client.get('/nonexistent')
        
        assert response.status_code == 404
        
        data = json.loads(response.data)
        assert data['success'] is False
        assert 'Endpoint not found' in data['error']
    
    @patch('app.pubchem_client.search_compound')
    def test_internal_server_error(self, mock_search, client):
        """Test 500 error handling."""
        mock_search.side_effect = Exception("Unexpected error")
        
        response = client.post('/api/pubchem/search',
                             json={'query': 'test'},
                             content_type='application/json')
        
        assert response.status_code == 500
        
        data = json.loads(response.data)
        assert data['success'] is False
        assert 'Internal server error' in data['error']