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
    get_supported_parameters
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


if __name__ == "__main__":
    pytest.main([__file__])