#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
PySCF Native App - 環境検証スクリプト
このスクリプトは、開発環境の健全性を検証します
"""

import sys
import importlib
import subprocess
from pathlib import Path
from typing import List, Tuple, Dict, Any

# カラー出力用の定数
class Colors:
    RED = '\033[0;31m'
    GREEN = '\033[0;32m'
    YELLOW = '\033[1;33m'
    BLUE = '\033[0;34m'
    PURPLE = '\033[0;35m'
    CYAN = '\033[0;36m'
    NC = '\033[0m'  # No Color

def log_info(message: str) -> None:
    """情報ログを出力"""
    print(f"{Colors.BLUE}[INFO]{Colors.NC} {message}")

def log_success(message: str) -> None:
    """成功ログを出力"""
    print(f"{Colors.GREEN}[SUCCESS]{Colors.NC} {message}")

def log_warning(message: str) -> None:
    """警告ログを出力"""
    print(f"{Colors.YELLOW}[WARNING]{Colors.NC} {message}")

def log_error(message: str) -> None:
    """エラーログを出力"""
    print(f"{Colors.RED}[ERROR]{Colors.NC} {message}")

def check_python_version() -> bool:
    """Python バージョンをチェック"""
    log_info(f"Python バージョンをチェック中...")
    version = sys.version_info
    version_str = f"{version.major}.{version.minor}.{version.micro}"
    
    if version.major == 3 and version.minor >= 10:
        log_success(f"Python {version_str} ✓")
        return True
    else:
        log_error(f"Python {version_str} - Python 3.10+ が必要です")
        return False

def check_required_packages() -> bool:
    """必須パッケージの存在確認"""
    required_packages = [
        ('pyscf', 'PySCF - 量子化学計算'),
        ('rdkit', 'RDKit - 化学情報学'),
        ('geometric', 'geometric - 分子幾何最適化'),
        ('flask', 'Flask - ウェブフレームワーク'),
        ('flask_cors', 'Flask-CORS - CORS対応'),
        ('pydantic', 'Pydantic - データバリデーション'),
        ('gevent', 'Gevent - 非同期処理'),
        ('requests', 'Requests - HTTP クライアント'),
        ('conda_pack', 'conda-pack - 環境パッケージ化'),
        ('datamodel_code_generator', 'datamodel-code-generator - コード生成'),
        ('PyInstaller', 'PyInstaller - Python実行ファイル作成'),
    ]
    
    log_info("必須パッケージをチェック中...")
    all_available = True
    
    for package_name, description in required_packages:
        try:
            module = importlib.import_module(package_name)
            version = getattr(module, '__version__', 'バージョン不明')
            log_success(f"{description}: {version} ✓")
        except ImportError:
            log_error(f"{description}: インストールされていません ✗")
            all_available = False
    
    return all_available

def check_pyscf_functionality() -> bool:
    """PySCF の基本機能をテスト"""
    log_info("PySCF の基本機能をテスト中...")
    
    try:
        import pyscf
        from pyscf import gto, scf
        
        # 水分子の簡単な計算テスト
        mol = gto.Mole()
        mol.atom = 'H 0 0 0; H 0 0 1.1'
        mol.basis = 'sto-3g'
        mol.build()
        
        mf = scf.RHF(mol)
        energy = mf.kernel()
        
        if abs(energy - (-1.0671)) < 0.1:  # 大まかな期待値チェック
            log_success("PySCF 計算テスト成功 ✓")
            return True
        else:
            log_warning(f"PySCF 計算結果が予期と異なります: {energy}")
            return False
            
    except Exception as e:
        log_error(f"PySCF テストに失敗: {e}")
        return False

def check_rdkit_functionality() -> bool:
    """RDKit の基本機能をテスト"""
    log_info("RDKit の基本機能をテスト中...")
    
    try:
        from rdkit import Chem
        from rdkit.Chem import rdMolDescriptors
        
        # SMILES文字列からの分子生成テスト
        mol = Chem.MolFromSmiles('CCO')  # エタノール
        if mol is None:
            log_error("SMILES からの分子生成に失敗")
            return False
        
        # 分子量計算テスト
        mw = rdMolDescriptors.CalcExactMolWt(mol)
        expected_mw = 46.0419  # エタノールの分子量
        
        if abs(mw - expected_mw) < 0.01:
            log_success("RDKit 機能テスト成功 ✓")
            return True
        else:
            log_warning(f"RDKit 計算結果が予期と異なります: {mw}")
            return False
            
    except Exception as e:
        log_error(f"RDKit テストに失敗: {e}")
        return False

def check_flask_functionality() -> bool:
    """Flask の基本機能をテスト"""
    log_info("Flask の基本機能をテスト中...")
    
    try:
        from flask import Flask
        from flask_cors import CORS
        import flask_sock
        
        # 簡単なFlaskアプリケーション作成テスト
        app = Flask(__name__)
        CORS(app)
        
        @app.route('/test')
        def test():
            return {'message': 'Flask test successful'}
        
        # テストクライアントで簡単なテスト
        with app.test_client() as client:
            response = client.get('/test')
            if response.status_code == 200:
                log_success("Flask 機能テスト成功 ✓")
                return True
            else:
                log_error(f"Flask テスト失敗: ステータスコード {response.status_code}")
                return False
                
    except Exception as e:
        log_error(f"Flask テストに失敗: {e}")
        return False

def check_conda_environment() -> bool:
    """conda 環境の確認"""
    log_info("conda 環境をチェック中...")
    
    try:
        # conda 環境名の確認
        conda_env = subprocess.run(
            ['conda', 'info', '--json'],
            capture_output=True,
            text=True,
            check=True
        )
        
        import json
        env_info = json.loads(conda_env.stdout)
        active_env = env_info.get('active_prefix_name', 'base')
        
        if active_env == 'pyscf-env':
            log_success(f"適切なconda環境がアクティブです: {active_env} ✓")
            return True
        else:
            log_warning(f"現在のconda環境: {active_env} (pyscf-env が推奨)")
            return False
            
    except subprocess.CalledProcessError:
        log_warning("conda コマンドが利用できません")
        return False
    except Exception as e:
        log_error(f"conda 環境チェック失敗: {e}")
        return False

def check_project_structure() -> bool:
    """プロジェクト構造の確認"""
    log_info("プロジェクト構造をチェック中...")
    
    script_dir = Path(__file__).parent
    project_root = script_dir.parent
    
    required_files = [
        'package.json',
        'tsconfig.json',
        'webpack.config.ts',
        '.github/environment.yml',
        'src/main.ts',
        'src/python/app.py',
        'src/api-spec/openapi.yaml'
    ]
    
    all_exist = True
    for file_path in required_files:
        full_path = project_root / file_path
        if full_path.exists():
            log_success(f"{file_path} ✓")
        else:
            log_error(f"{file_path} が見つかりません ✗")
            all_exist = False
    
    return all_exist

def main() -> None:
    """メイン実行関数"""
    print(f"{Colors.CYAN}=== PySCF Native App 環境検証 ==={Colors.NC}")
    print()
    
    tests = [
        ("Python バージョン", check_python_version),
        ("必須パッケージ", check_required_packages),
        ("PySCF 機能", check_pyscf_functionality),
        ("RDKit 機能", check_rdkit_functionality),
        ("Flask 機能", check_flask_functionality),
        ("conda 環境", check_conda_environment),
        ("プロジェクト構造", check_project_structure),
    ]
    
    results = []
    for test_name, test_func in tests:
        print(f"\n{Colors.PURPLE}--- {test_name} ---{Colors.NC}")
        try:
            result = test_func()
            results.append((test_name, result))
        except Exception as e:
            log_error(f"{test_name}テスト中にエラー: {e}")
            results.append((test_name, False))
    
    # 結果サマリー
    print(f"\n{Colors.CYAN}=== 検証結果サマリー ==={Colors.NC}")
    passed = sum(1 for _, result in results if result)
    total = len(results)
    
    for test_name, result in results:
        status = f"{Colors.GREEN}PASS{Colors.NC}" if result else f"{Colors.RED}FAIL{Colors.NC}"
        print(f"  {test_name}: {status}")
    
    print(f"\n総合結果: {passed}/{total} テスト合格")
    
    if passed == total:
        log_success("すべての検証が成功しました！環境は正常です。")
        sys.exit(0)
    else:
        log_error("いくつかの検証が失敗しました。環境を確認してください。")
        sys.exit(1)

if __name__ == '__main__':
    main()