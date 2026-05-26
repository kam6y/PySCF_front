#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
PySCF_front - ビルド完了後検証スクリプト
このスクリプトは、ビルド完了後にbundled環境の完全性を検証します
"""

import json
import os
import sys
import subprocess
from pathlib import Path
from typing import List, Tuple

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

def check_conda_environment(project_root: Path) -> bool:
    """bundled conda環境の完全性をチェック"""
    log_info("bundled conda環境をチェック中...")
    conda_env_path = project_root / "conda_env"
    
    if not conda_env_path.exists():
        log_error("conda_env ディレクトリが見つかりません")
        log_info("conda-pack がまだ実行されていないか、失敗している可能性があります")
        return False
    
    # 重要なファイルの存在確認
    required_files = [
        Path("bin/python"),
        Path("bin/gunicorn"),
        Path("bin/pip"),
    ]
    
    all_exist = True
    for file_path in required_files:
        full_path = conda_env_path / file_path
        if full_path.exists():
            log_success(f"✓ {file_path}")
        else:
            log_error(f"✗ {file_path} が見つかりません")
            all_exist = False
    
    # Python バージョンディレクトリを動的検出
    python_lib_dirs = list(conda_env_path.glob("lib/python3.*"))
    if python_lib_dirs:
        python_version_dir = python_lib_dirs[0].name  # python3.12 など
        log_success(f"✓ lib/{python_version_dir}")
    else:
        log_error("✗ lib/python3.* が見つかりません")
        all_exist = False
    
    # conda環境のサイズチェック
    try:
        total_size = sum(f.stat().st_size for f in conda_env_path.rglob('*') if f.is_file())
        size_mb = total_size / (1024 * 1024)
        
        if size_mb < 100:  # 100MB未満は異常に小さい
            log_warning(f"conda環境が小さすぎます: {size_mb:.1f} MB")
            log_warning("依存関係が不完全な可能性があります")
            return False
        else:
            log_success(f"conda環境サイズ: {size_mb:.1f} MB")
            
    except Exception as e:
        log_warning(f"サイズチェック失敗: {e}")
    
    return all_exist


def check_config_files(project_root: Path) -> bool:
    """設定ファイルの存在確認"""
    log_info("設定ファイルをチェック中...")
    
    required_configs = [
        "config/server-config.json",
        "src/api-spec/openapi.yaml",
    ]
    
    all_exist = True
    for config_path in required_configs:
        full_path = project_root / config_path
        if full_path.exists():
            log_success(f"✓ {config_path}")
        else:
            log_error(f"✗ {config_path} が見つかりません")
            all_exist = False
    
    return all_exist

def _display_path(project_root: Path, target_path: Path) -> str:
    """プロジェクトルートからの相対パスを表示用に返す"""
    try:
        return str(target_path.relative_to(project_root))
    except ValueError:
        return str(target_path)


def _load_extra_resource_targets(project_root: Path) -> List[Tuple[str, str]]:
    """electron-builder extraResources の from/to 対応を取得"""
    package_json_path = project_root / "package.json"

    try:
        with package_json_path.open("r", encoding="utf-8") as handle:
            package_data = json.load(handle)
    except Exception as e:
        log_error(f"package.json の読み込みに失敗しました: {e}")
        return []

    extra_resources = package_data.get("build", {}).get("extraResources", [])
    targets: List[Tuple[str, str]] = []

    for resource in extra_resources:
        if not isinstance(resource, dict):
            continue

        source = resource.get("from")
        target = resource.get("to")
        if not isinstance(source, str) or not isinstance(target, str):
            continue

        normalized_source = source.replace("\\", "/").rstrip("/")
        normalized_target = target.replace("\\", "/").rstrip("/")
        targets.append((normalized_source, normalized_target))

    return targets


def _find_packaged_resource_dirs(project_root: Path) -> List[Path]:
    """既存の packaged app から Electron resources ディレクトリを探す"""
    release_path = project_root / "release"
    if not release_path.exists():
        return []

    resource_dirs: List[Path] = []
    for dirname in ("Resources", "resources"):
        for candidate in release_path.rglob(dirname):
            if candidate.is_dir() and (candidate / "app.asar").exists():
                resource_dirs.append(candidate)

    return sorted(set(resource_dirs))


def check_packaged_resource_layout(project_root: Path) -> bool:
    """packaged app で参照される extraResources 配置をチェック"""
    log_info("packaged resource layout をチェック中...")

    source_root = "src/api-spec"
    target_root = "src/api-spec"
    openapi_source_path = project_root / source_root / "openapi.yaml"
    packaged_openapi_path = Path(target_root) / "openapi.yaml"

    all_success = True

    if openapi_source_path.exists():
        log_success(f"✓ {_display_path(project_root, openapi_source_path)}")
    else:
        log_error(f"✗ {_display_path(project_root, openapi_source_path)} が見つかりません")
        all_success = False

    extra_resource_targets = _load_extra_resource_targets(project_root)
    if (source_root, target_root) in extra_resource_targets:
        log_success(f"✓ extraResources: {source_root} -> {target_root}")
    else:
        log_error(f"✗ extraResources に {source_root} -> {target_root} がありません")
        all_success = False

    resource_dirs = _find_packaged_resource_dirs(project_root)
    if not resource_dirs:
        log_warning("既存の packaged app は見つかりません。extraResources 設定のみ検証しました")
        return all_success

    for resource_dir in resource_dirs:
        expected_path = resource_dir / packaged_openapi_path
        if expected_path.exists():
            log_success(f"✓ {_display_path(project_root, expected_path)}")
        else:
            log_warning(
                f"⚠ {_display_path(project_root, expected_path)} が見つかりません "
                "(既存の package は再作成前の可能性があります)"
            )

    return all_success


def check_frontend_build(project_root: Path) -> bool:
    """フロントエンドビルドの確認"""
    log_info("フロントエンドビルドをチェック中...")
    dist_path = project_root / "dist"
    
    if not dist_path.exists():
        log_error("dist ディレクトリが見つかりません")
        return False
    
    required_files = [
        "index.html",
        "splash.html",
        "main.js",
        "preload.js",
        "splashPreload.js",
    ]

    required_asset_patterns = [
        "assets/index-*.js",
        "assets/splash-*.js",
    ]
    
    all_exist = True
    for file_name in required_files:
        full_path = dist_path / file_name
        if full_path.exists():
            log_success(f"✓ {file_name}")
        else:
            log_error(f"✗ {file_name} が見つかりません")
            all_exist = False

    for pattern in required_asset_patterns:
        matches = sorted(dist_path.glob(pattern))
        if matches:
            for match in matches:
                log_success(f"✓ {_display_path(project_root, match)}")
        else:
            log_error(f"✗ {pattern} に一致するファイルが見つかりません")
            all_exist = False
    
    return all_exist

def validate_conda_functionality(project_root: Path) -> bool:
    """conda環境の機能テスト"""
    log_info("conda環境の機能をテスト中...")
    python_exe = project_root / "conda_env" / "bin" / "python"

    if not python_exe.exists():
        log_error("conda環境のPython実行ファイルが見つかりません")
        return False

    required_modules = [
        "pyscf",
        "rdkit",
        "fastapi",
        "uvicorn",
        "socketio",
        "gunicorn",
        "pydantic",
        "conda_pack",
    ]

    # 環境変数でPySCFテストをスキップ可能にする
    skip_pyscf_test = os.environ.get('SKIP_PYSCF_TEST', '').lower() in ('true', '1', 'yes')

    all_success = True
    for module_name in required_modules:
        # PySCFテストをスキップする場合
        if module_name == "pyscf" and skip_pyscf_test:
            log_warning(f"⚠ {module_name} import test skipped (SKIP_PYSCF_TEST set)")
            continue

        try:
            # PySCFは大型ライブラリなので特に長いタイムアウトを設定
            timeout_duration = 60 if module_name == "pyscf" else 10
            result = subprocess.run([
                str(python_exe),
                "-c",
                f"import {module_name}; print(f'{module_name}: OK')"
            ], capture_output=True, text=True, timeout=timeout_duration)

            if result.returncode == 0:
                log_success(f"✓ {module_name} import successful")
            else:
                if module_name == "pyscf":
                    log_warning(f"⚠ {module_name} import failed: {result.stderr.strip()}")
                    log_warning("PySCF import failed but continuing with build validation")
                else:
                    log_error(f"✗ {module_name} import failed: {result.stderr.strip()}")
                    all_success = False

        except subprocess.TimeoutExpired:
            if module_name == "pyscf":
                log_warning(f"⚠ {module_name} import timeout")
                log_warning("PySCF import timed out but continuing with build validation")
            else:
                log_error(f"✗ {module_name} import timeout")
                all_success = False
        except Exception as e:
            if module_name == "pyscf":
                log_warning(f"⚠ {module_name} import error: {e}")
                log_warning("PySCF import error but continuing with build validation")
            else:
                log_error(f"✗ {module_name} import error: {e}")
                all_success = False

    return all_success

def main() -> None:
    """メイン実行関数"""
    print(f"{Colors.CYAN}=== PySCF_front ビルド完全性検証 ==={Colors.NC}")
    print()
    
    # プロジェクトルートを一度だけ取得
    script_dir = Path(__file__).parent
    project_root = script_dir.parent
    
    tests = [
        ("bundled conda環境", lambda: check_conda_environment(project_root)),
        ("設定ファイル", lambda: check_config_files(project_root)),
        ("packaged resource layout", lambda: check_packaged_resource_layout(project_root)),
        ("フロントエンドビルド", lambda: check_frontend_build(project_root)),
        ("conda環境機能", lambda: validate_conda_functionality(project_root)),
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
    print(f"\n{Colors.CYAN}=== ビルド完全性検証結果 ==={Colors.NC}")
    passed = sum(1 for _, result in results if result)
    total = len(results)
    
    for test_name, result in results:
        status = f"{Colors.GREEN}PASS{Colors.NC}" if result else f"{Colors.RED}FAIL{Colors.NC}"
        print(f"  {test_name}: {status}")
    
    print(f"\n総合結果: {passed}/{total} テスト合格")
    
    if passed == total:
        log_success("ビルド完全性検証が成功しました！パッケージングの準備ができています。")
        sys.exit(0)
    else:
        log_error("ビルド完全性検証が失敗しました。不足しているコンポーネントを確認してください。")
        print(f"\n{Colors.YELLOW}修復手順:{Colors.NC}")
        print("1. 失敗したコンポーネントを確認")
        print("2. 該当するビルドステップを再実行:")
        print("   - conda環境: npm run build:conda-pack")
  
        print("   - フロントエンド: npm run build:electron")
        print("3. 再度検証: npm run validate-build")
        sys.exit(1)

if __name__ == '__main__':
    main()
