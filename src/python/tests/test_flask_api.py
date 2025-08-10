"""Tests for Flask API endpoints."""

import pytest
import json
from unittest.mock import Mock, patch

# Add the parent directory to the path for imports
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from app import app
from pubchem.client import PubChemError, CompoundData

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
    
    #
    # 変更点: モックを簡素化し、テストの保守性を向上
    #
    @patch('app.pubchem_client.search_compound')
    def test_search_pubchem_success(self, mock_search_compound, client):
        """Test successful PubChem search."""
        # Mock the entire CompoundData object that the service method returns
        mock_compound_data = CompoundData(
            cid=962,
            iupac_name='oxidane',
            molecular_formula='H2O',
            molecular_weight=18.015,
            synonyms=['Water', 'H2O'],
            atoms=[
                {'element': 'O', 'x': 0.0, 'y': 0.0, 'z': 0.1193},
                {'element': 'H', 'x': 0.0, 'y': 0.7632, 'z': -0.477},
                {'element': 'H', 'x': 0.0, 'y': -0.7632, 'z': -0.477}
            ]
        )
        mock_search_compound.return_value = mock_compound_data
        
        # Make request
        response = client.post('/api/pubchem/search',
                             json={'query': 'water', 'search_type': 'name'},
                             content_type='application/json')
        
        assert response.status_code == 200
        
        data = json.loads(response.data)
        assert data['success'] is True
        assert 'xyz' in data['data']
        assert '3\noxidane (H2O) - CID: 962' in data['data']['xyz']
        assert data['data']['compound_info']['cid'] == 962
        assert data['data']['compound_info']['iupac_name'] == 'oxidane'
        assert data['data']['atom_count'] == 3

    # (以降のテストは変更なし)
    # ...