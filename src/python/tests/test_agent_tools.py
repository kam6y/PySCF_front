"""
AIエージェント用ツールラッパーのテスト

agent/tools.pyモジュールの各ツール関数のユニットテストを実装。
正常系、異常系、エラーハンドリングを包括的にテストします。
"""

import json
import pytest
from unittest.mock import Mock, patch, MagicMock

import requests

# テスト対象のモジュールをインポート
from agent.tools import (
    list_all_calculations,
    get_calculation_details,
    get_supported_parameters,
    start_quantum_calculation,
    search_pubchem_by_name
)


class TestAgentTools:
    """AIエージェント用ツールラッパーのテストクラス"""
    
    def test_list_all_calculations_success(self):
        """list_all_calculations関数の正常系テスト"""
        # モックレスポンスデータ
        mock_response_data = {
            "success": True,
            "data": {
                "base_directory": "/path/to/calculations",
                "calculations": [
                    {
                        "id": "calc_20240101_120000_abcd1234",
                        "name": "Water DFT Calculation",
                        "path": "/path/to/calc1",
                        "date": "2024-01-01 12:00:00",
                        "has_checkpoint": True,
                        "status": "completed"
                    },
                    {
                        "id": "calc_20240102_130000_efgh5678",
                        "name": "Benzene Optimization",
                        "path": "/path/to/calc2",
                        "date": "2024-01-02 13:00:00",
                        "has_checkpoint": True,
                        "status": "running"
                    }
                ],
                "count": 2
            }
        }
        
        # requests.getをモック
        with patch('agent.tools.requests.get') as mock_get:
            mock_response = Mock()
            mock_response.json.return_value = mock_response_data
            mock_response.raise_for_status.return_value = None
            mock_get.return_value = mock_response
            
            # 関数実行
            result = list_all_calculations()
            
            # 結果検証
            assert isinstance(result, str)
            parsed_result = json.loads(result)
            assert parsed_result["success"] is True
            assert parsed_result["data"]["count"] == 2
            assert len(parsed_result["data"]["calculations"]) == 2
            
            # API呼び出し確認
            mock_get.assert_called_once_with(
                "http://127.0.0.1:5000/api/quantum/calculations",
                timeout=30
            )
    
    def test_list_all_calculations_http_error(self):
        """list_all_calculations関数のHTTPエラーテスト"""
        with patch('agent.tools.requests.get') as mock_get:
            # 404エラーをシミュレート
            mock_response = Mock()
            mock_response.status_code = 404
            mock_get.return_value = mock_response
            mock_get.return_value.raise_for_status.side_effect = requests.exceptions.HTTPError(response=mock_response)
            
            result = list_all_calculations()
            
            assert "Error:" in result
            assert "計算リストが見つかりませんでした" in result
    
    def test_list_all_calculations_connection_error(self):
        """list_all_calculations関数の接続エラーテスト"""
        with patch('agent.tools.requests.get') as mock_get:
            mock_get.side_effect = requests.exceptions.ConnectionError("Connection failed")
            
            result = list_all_calculations()
            
            assert "Error:" in result
            assert "APIサーバーに接続できませんでした" in result
    
    def test_list_all_calculations_timeout_error(self):
        """list_all_calculations関数のタイムアウトエラーテスト"""
        with patch('agent.tools.requests.get') as mock_get:
            mock_get.side_effect = requests.exceptions.Timeout("Request timed out")
            
            result = list_all_calculations()
            
            assert "Error:" in result
            assert "タイムアウトしました" in result
    
    def test_get_calculation_details_success(self):
        """get_calculation_details関数の正常系テスト"""
        test_calc_id = "calc_20240101_120000_abcd1234"
        mock_response_data = {
            "success": True,
            "data": {
                "calculation": {
                    "id": test_calc_id,
                    "name": "Water DFT Calculation",
                    "status": "completed",
                    "parameters": {
                        "calculation_method": "DFT",
                        "basis_function": "6-31G(d)",
                        "exchange_correlation": "B3LYP",
                        "charges": 0,
                        "spin": 0
                    },
                    "results": {
                        "scf_energy": -76.3680423,
                        "converged": True,
                        "homo_index": 4,
                        "lumo_index": 5
                    }
                }
            }
        }
        
        with patch('agent.tools.requests.get') as mock_get:
            mock_response = Mock()
            mock_response.json.return_value = mock_response_data
            mock_response.raise_for_status.return_value = None
            mock_get.return_value = mock_response
            
            result = get_calculation_details(test_calc_id)
            
            assert isinstance(result, str)
            parsed_result = json.loads(result)
            assert parsed_result["success"] is True
            assert parsed_result["data"]["calculation"]["id"] == test_calc_id
            assert parsed_result["data"]["calculation"]["results"]["scf_energy"] == -76.3680423
            
            # API呼び出し確認
            mock_get.assert_called_once_with(
                f"http://127.0.0.1:5000/api/quantum/calculations/{test_calc_id}",
                timeout=30
            )
    
    def test_get_calculation_details_invalid_id(self):
        """get_calculation_details関数の無効なIDテスト"""
        # 空文字列の場合
        result = get_calculation_details("")
        assert "Error:" in result
        assert "calculation_idは必須の文字列パラメータです" in result
        
        # Noneの場合
        result = get_calculation_details(None)
        assert "Error:" in result
        assert "calculation_idは必須の文字列パラメータです" in result
    
    def test_get_calculation_details_not_found(self):
        """get_calculation_details関数の計算が見つからないテスト"""
        test_calc_id = "invalid_calc_id"
        
        with patch('agent.tools.requests.get') as mock_get:
            mock_response = Mock()
            mock_response.status_code = 404
            mock_get.return_value = mock_response
            mock_get.return_value.raise_for_status.side_effect = requests.exceptions.HTTPError(response=mock_response)
            
            result = get_calculation_details(test_calc_id)
            
            assert "Error:" in result
            assert f"ID '{test_calc_id}' の計算が見つかりませんでした" in result
    
    def test_get_calculation_details_server_error(self):
        """get_calculation_details関数のサーバーエラーテスト"""
        test_calc_id = "calc_20240101_120000_abcd1234"
        
        with patch('agent.tools.requests.get') as mock_get:
            mock_response = Mock()
            mock_response.status_code = 500
            mock_get.return_value = mock_response
            mock_get.return_value.raise_for_status.side_effect = requests.exceptions.HTTPError(response=mock_response)
            
            result = get_calculation_details(test_calc_id)
            
            assert "Error:" in result
            assert "サーバー内部エラーが発生しました" in result
    
    def test_get_supported_parameters_success(self):
        """get_supported_parameters関数の正常系テスト"""
        mock_response_data = {
            "success": True,
            "data": {
                "calculation_methods": ["DFT", "HF", "MP2", "CCSD", "TDDFT", "CASCI", "CASSCF"],
                "basis_functions": {
                    "Minimal": ["STO-3G"],
                    "Pople Style": ["6-31G", "6-31G(d)", "6-31+G(d,p)", "6-311G(d,p)"],
                    "Correlation Consistent": ["cc-pVDZ", "cc-pVTZ", "aug-cc-pVDZ", "aug-cc-pVTZ"]
                },
                "exchange_correlation": {
                    "Hybrid": ["B3LYP", "PBE0", "M06-2X", "CAM-B3LYP"],
                    "GGA": ["PBE", "BLYP"],
                    "Meta-GGA": ["M06", "TPSS"]
                },
                "solvent_methods": ["none", "ief-pcm", "c-pcm", "cosmo"],
                "solvents": {
                    "Highly Polar": [
                        {"value": "water", "display": "Water", "dielectric_constant": 78.39},
                        {"value": "dimethylsulfoxide", "display": "DMSO", "dielectric_constant": 46.83}
                    ]
                },
                "tddft_methods": ["TDDFT", "TDA"]
            }
        }
        
        with patch('agent.tools.requests.get') as mock_get:
            mock_response = Mock()
            mock_response.json.return_value = mock_response_data
            mock_response.raise_for_status.return_value = None
            mock_get.return_value = mock_response
            
            result = get_supported_parameters()
            
            assert isinstance(result, str)
            parsed_result = json.loads(result)
            assert parsed_result["success"] is True
            assert "DFT" in parsed_result["data"]["calculation_methods"]
            assert "B3LYP" in parsed_result["data"]["exchange_correlation"]["Hybrid"]
            
            # API呼び出し確認
            mock_get.assert_called_once_with(
                "http://127.0.0.1:5000/api/quantum/supported-parameters",
                timeout=30
            )
    
    def test_get_supported_parameters_server_error(self):
        """get_supported_parameters関数のサーバーエラーテスト"""
        with patch('agent.tools.requests.get') as mock_get:
            mock_response = Mock()
            mock_response.status_code = 500
            mock_get.return_value = mock_response
            mock_get.return_value.raise_for_status.side_effect = requests.exceptions.HTTPError(response=mock_response)
            
            result = get_supported_parameters()
            
            assert "Error:" in result
            assert "サーバー内部エラーが発生しました" in result
    
    def test_get_supported_parameters_json_decode_error(self):
        """get_supported_parameters関数のJSON解析エラーテスト"""
        with patch('agent.tools.requests.get') as mock_get:
            mock_response = Mock()
            mock_response.raise_for_status.return_value = None
            mock_response.json.side_effect = json.JSONDecodeError("Invalid JSON", "doc", 0)
            mock_get.return_value = mock_response
            
            result = get_supported_parameters()
            
            assert "Error:" in result
            assert "APIレスポンスの解析に失敗しました" in result
    
    def test_all_functions_timeout_handling(self):
        """全ての関数のタイムアウト処理テスト"""
        functions_to_test = [
            (list_all_calculations, []),
            (get_calculation_details, ["test_id"]),
            (get_supported_parameters, [])
        ]
        
        for func, args in functions_to_test:
            with patch('agent.tools.requests.get') as mock_get:
                mock_get.side_effect = requests.exceptions.Timeout("Request timed out")
                
                result = func(*args)
                
                assert "Error:" in result
                assert "タイムアウトしました" in result
    
    def test_all_functions_unexpected_error_handling(self):
        """全ての関数の予期しないエラー処理テスト"""
        functions_to_test = [
            (list_all_calculations, []),
            (get_calculation_details, ["test_id"]),
            (get_supported_parameters, [])
        ]
        
        for func, args in functions_to_test:
            with patch('agent.tools.requests.get') as mock_get:
                mock_get.side_effect = Exception("Unexpected error")
                
                result = func(*args)
                
                assert "Error:" in result
                assert "予期しないエラーが発生しました" in result
    
    def test_start_quantum_calculation_success(self):
        """start_quantum_calculation関数の正常系テスト"""
        test_xyz = "3\nWater molecule\nO 0.0 0.0 0.0\nH 0.7570 0.5860 0.0\nH -0.7570 0.5860 0.0"
        mock_response_data = {
            "success": True,
            "data": {
                "calculation": {
                    "id": "calc_20240101_120000_abcd1234",
                    "name": "AI Generated Calculation",
                    "status": "running",
                    "parameters": {
                        "calculation_method": "DFT",
                        "basis_function": "6-31G(d)",
                        "exchange_correlation": "B3LYP",
                        "charges": 0,
                        "spin": 0
                    },
                    "created_at": "2024-01-01T12:00:00Z"
                }
            }
        }
        
        with patch('agent.tools.requests.post') as mock_post:
            mock_response = Mock()
            mock_response.json.return_value = mock_response_data
            mock_response.raise_for_status.return_value = None
            mock_post.return_value = mock_response
            
            result = start_quantum_calculation(
                xyz=test_xyz,
                calculation_method="DFT",
                basis_function="6-31G(d)",
                charges=0,
                spin=0,
                exchange_correlation="B3LYP",
                name="Test Calculation"
            )
            
            assert isinstance(result, str)
            parsed_result = json.loads(result)
            assert parsed_result["success"] is True
            assert parsed_result["data"]["calculation"]["status"] == "running"
            
            # API呼び出し確認
            mock_post.assert_called_once()
            call_args = mock_post.call_args
            assert call_args[0][0] == "http://127.0.0.1:5000/api/quantum/calculate"
            assert call_args[1]["json"]["xyz"] == test_xyz.strip()
            assert call_args[1]["json"]["calculation_method"] == "DFT"
    
    def test_start_quantum_calculation_validation_errors(self):
        """start_quantum_calculation関数のバリデーションエラーテスト"""
        # XYZパラメータが空の場合
        result = start_quantum_calculation(xyz="")
        assert "Error:" in result
        assert "xyz parameter is required" in result
        
        # XYZパラメータがNoneの場合
        result = start_quantum_calculation(xyz=None)
        assert "Error:" in result
        assert "xyz parameter is required" in result
        
        # 無効な計算メソッドの場合
        result = start_quantum_calculation(xyz="test", calculation_method="INVALID")
        assert "Error:" in result
        assert "calculation_method must be one of" in result
        
        # 無効な電荷の場合
        result = start_quantum_calculation(xyz="test", charges=-20)
        assert "Error:" in result
        assert "charges must be an integer between -10 and 10" in result
        
        # 無効なスピンの場合
        result = start_quantum_calculation(xyz="test", spin=-1)
        assert "Error:" in result
        assert "spin must be an integer between 0 and 10" in result
        
        # 無効な名前の場合
        result = start_quantum_calculation(xyz="test", name="")
        assert "Error:" in result
        assert "name must be a non-empty string" in result
    
    def test_start_quantum_calculation_http_errors(self):
        """start_quantum_calculation関数のHTTPエラーテスト"""
        test_xyz = "3\nWater\nO 0.0 0.0 0.0\nH 1.0 0.0 0.0\nH 0.0 1.0 0.0"
        
        # 400エラー（バリデーションエラー）
        with patch('agent.tools.requests.post') as mock_post:
            mock_response = Mock()
            mock_response.status_code = 400
            mock_response.json.return_value = {"message": "Invalid XYZ format"}
            mock_post.return_value = mock_response
            mock_post.return_value.raise_for_status.side_effect = requests.exceptions.HTTPError(response=mock_response)
            
            result = start_quantum_calculation(xyz=test_xyz)
            
            assert "Error:" in result
            assert "Invalid calculation parameters" in result
        
        # 500エラー（サーバーエラー）
        with patch('agent.tools.requests.post') as mock_post:
            mock_response = Mock()
            mock_response.status_code = 500
            mock_post.return_value = mock_response
            mock_post.return_value.raise_for_status.side_effect = requests.exceptions.HTTPError(response=mock_response)
            
            result = start_quantum_calculation(xyz=test_xyz)
            
            assert "Error:" in result
            assert "Internal server error occurred while starting calculation" in result
    
    def test_start_quantum_calculation_network_errors(self):
        """start_quantum_calculation関数のネットワークエラーテスト"""
        test_xyz = "3\nWater\nO 0.0 0.0 0.0\nH 1.0 0.0 0.0\nH 0.0 1.0 0.0"
        
        # 接続エラー
        with patch('agent.tools.requests.post') as mock_post:
            mock_post.side_effect = requests.exceptions.ConnectionError("Connection failed")
            
            result = start_quantum_calculation(xyz=test_xyz)
            
            assert "Error:" in result
            assert "Could not connect to API server" in result
        
        # タイムアウトエラー
        with patch('agent.tools.requests.post') as mock_post:
            mock_post.side_effect = requests.exceptions.Timeout("Request timed out")
            
            result = start_quantum_calculation(xyz=test_xyz)
            
            assert "Error:" in result
            assert "API request timed out" in result
    
    def test_search_pubchem_by_name_success(self):
        """search_pubchem_by_name関数の正常系テスト"""
        test_compound = "water"
        mock_response_data = {
            "success": True,
            "data": {
                "compound": {
                    "name": "Water",
                    "cid": 962,
                    "molecular_formula": "H2O",
                    "molecular_weight": 18.015,
                    "xyz_structure": "3\nWater\nO 0.000000 0.000000 0.119159\nH 0.000000 0.757000 -0.476637\nH 0.000000 -0.757000 -0.476637",
                    "structure_source": "PubChem 3D conformer"
                }
            }
        }
        
        with patch('agent.tools.requests.post') as mock_post:
            mock_response = Mock()
            mock_response.json.return_value = mock_response_data
            mock_response.raise_for_status.return_value = None
            mock_post.return_value = mock_response
            
            result = search_pubchem_by_name(compound_name=test_compound)
            
            assert isinstance(result, str)
            parsed_result = json.loads(result)
            assert parsed_result["success"] is True
            assert parsed_result["data"]["compound"]["name"] == "Water"
            assert parsed_result["data"]["compound"]["cid"] == 962
            
            # API呼び出し確認
            mock_post.assert_called_once()
            call_args = mock_post.call_args
            assert call_args[0][0] == "http://127.0.0.1:5000/api/pubchem/search"
            assert call_args[1]["json"]["query"] == test_compound
            assert call_args[1]["json"]["searchType"] == "name"
    
    def test_search_pubchem_by_name_validation_errors(self):
        """search_pubchem_by_name関数のバリデーションエラーテスト"""
        # 空のcompound_nameの場合
        result = search_pubchem_by_name(compound_name="")
        assert "Error:" in result
        assert "compound_name parameter is required" in result
        
        # Noneのcompound_nameの場合
        result = search_pubchem_by_name(compound_name=None)
        assert "Error:" in result
        assert "compound_name parameter is required" in result
        
        # 無効なsearch_typeの場合
        result = search_pubchem_by_name(compound_name="water", search_type="invalid")
        assert "Error:" in result
        assert "search_type must be one of" in result
        
        # 空白のみのcompound_nameの場合
        result = search_pubchem_by_name(compound_name="   ")
        assert "Error:" in result
        assert "compound_name cannot be empty or contain only whitespace" in result
    
    def test_search_pubchem_by_name_not_found(self):
        """search_pubchem_by_name関数の化合物が見つからないテスト"""
        test_compound = "nonexistentcompound12345"
        
        with patch('agent.tools.requests.post') as mock_post:
            mock_response = Mock()
            mock_response.status_code = 404
            mock_post.return_value = mock_response
            mock_post.return_value.raise_for_status.side_effect = requests.exceptions.HTTPError(response=mock_response)
            
            result = search_pubchem_by_name(compound_name=test_compound)
            
            assert "Error:" in result
            assert f"Compound '{test_compound}' not found in PubChem database" in result
    
    def test_search_pubchem_by_name_bad_request(self):
        """search_pubchem_by_name関数の不正なリクエストテスト"""
        test_compound = "invalid_search"
        
        with patch('agent.tools.requests.post') as mock_post:
            mock_response = Mock()
            mock_response.status_code = 400
            mock_response.json.return_value = {"message": "Invalid search query"}
            mock_post.return_value = mock_response
            mock_post.return_value.raise_for_status.side_effect = requests.exceptions.HTTPError(response=mock_response)
            
            result = search_pubchem_by_name(compound_name=test_compound)
            
            assert "Error:" in result
            assert "Invalid search request" in result
    
    def test_search_pubchem_by_name_different_search_types(self):
        """search_pubchem_by_name関数の異なる検索タイプテスト"""
        search_types = ["name", "cid", "formula"]
        
        for search_type in search_types:
            mock_response_data = {
                "success": True,
                "data": {"compound": {"name": "TestCompound"}}
            }
            
            with patch('agent.tools.requests.post') as mock_post:
                mock_response = Mock()
                mock_response.json.return_value = mock_response_data
                mock_response.raise_for_status.return_value = None
                mock_post.return_value = mock_response
                
                result = search_pubchem_by_name(
                    compound_name="test",
                    search_type=search_type
                )
                
                assert isinstance(result, str)
                parsed_result = json.loads(result)
                assert parsed_result["success"] is True
                
                # 正しい検索タイプが送信されているか確認
                call_args = mock_post.call_args
                assert call_args[1]["json"]["searchType"] == search_type
    
    def test_search_pubchem_by_name_network_errors(self):
        """search_pubchem_by_name関数のネットワークエラーテスト"""
        test_compound = "water"
        
        # 接続エラー
        with patch('agent.tools.requests.post') as mock_post:
            mock_post.side_effect = requests.exceptions.ConnectionError("Connection failed")
            
            result = search_pubchem_by_name(compound_name=test_compound)
            
            assert "Error:" in result
            assert "Could not connect to API server" in result
        
        # タイムアウトエラー
        with patch('agent.tools.requests.post') as mock_post:
            mock_post.side_effect = requests.exceptions.Timeout("Request timed out")
            
            result = search_pubchem_by_name(compound_name=test_compound)
            
            assert "Error:" in result
            assert "API request timed out" in result
        
        # JSON解析エラー
        with patch('agent.tools.requests.post') as mock_post:
            mock_response = Mock()
            mock_response.raise_for_status.return_value = None
            mock_response.json.side_effect = json.JSONDecodeError("Invalid JSON", "doc", 0)
            mock_post.return_value = mock_response
            
            result = search_pubchem_by_name(compound_name=test_compound)
            
            assert "Error:" in result
            assert "Failed to parse API response" in result
    
    def test_new_functions_integration_scenario(self):
        """新しい関数の統合シナリオテスト（PubChem検索→計算開始）"""
        # 1. PubChemで化合物を検索
        compound_search_response = {
            "success": True,
            "data": {
                "compound": {
                    "name": "Benzene",
                    "xyz_structure": "12\nBenzene\nC 1.4020 0.0000 0.0000\nC 0.7010 1.2135 0.0000\n..."
                }
            }
        }
        
        # 2. 計算開始
        calculation_start_response = {
            "success": True,
            "data": {
                "calculation": {
                    "id": "calc_20240101_120000_abcd1234",
                    "status": "running"
                }
            }
        }
        
        with patch('agent.tools.requests.post') as mock_post:
            # モックレスポンスを順番に設定
            mock_responses = [Mock(), Mock()]
            mock_responses[0].json.return_value = compound_search_response
            mock_responses[0].raise_for_status.return_value = None
            mock_responses[1].json.return_value = calculation_start_response
            mock_responses[1].raise_for_status.return_value = None
            
            mock_post.side_effect = mock_responses
            
            # 1. 化合物検索
            search_result = search_pubchem_by_name("benzene")
            search_data = json.loads(search_result)
            
            # 2. 検索結果から XYZ データを取得して計算開始
            xyz_data = search_data["data"]["compound"]["xyz_structure"]
            calc_result = start_quantum_calculation(
                xyz=xyz_data,
                calculation_method="DFT",
                name="Benzene DFT Calculation"
            )
            calc_data = json.loads(calc_result)
            
            # 結果検証
            assert search_data["success"] is True
            assert calc_data["success"] is True
            assert calc_data["data"]["calculation"]["status"] == "running"
            
            # 2回のAPI呼び出しが正しく行われたことを確認
            assert mock_post.call_count == 2


if __name__ == "__main__":
    pytest.main([__file__])