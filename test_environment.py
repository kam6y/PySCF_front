#!/usr/bin/env python
"""
開発環境テストスクリプト
"""
import sys
import os
from pathlib import Path

def test_python_version():
    """Python バージョンテスト"""
    print(f"Python version: {sys.version}")
    major, minor = sys.version_info[:2]
    if major >= 3 and minor >= 13:
        print("✓ Python version is compatible")
        return True
    else:
        print("✗ Python version should be 3.13 or higher")
        return False

def test_pyside6():
    """PySide6 インポートテスト"""
    try:
        from PySide6.QtWidgets import QApplication, QLabel
        from PySide6.QtCore import qVersion
        qt_version = qVersion()
        print(f"✓ PySide6 imported successfully (Qt {qt_version})")
        return True
    except ImportError as e:
        print(f"✗ PySide6 import failed: {e}")
        return False

def test_scientific_packages():
    """科学計算パッケージテスト"""
    packages = {
        'numpy': 'NumPy',
        'matplotlib': 'Matplotlib',
        'loguru': 'Loguru'
    }
    
    results = []
    for package_name, display_name in packages.items():
        try:
            __import__(package_name)
            print(f"✓ {display_name} imported successfully")
            results.append(True)
        except ImportError as e:
            print(f"✗ {display_name} import failed: {e}")
            results.append(False)
    
    return all(results)

def test_project_structure():
    """プロジェクト構造テスト"""
    required_dirs = [
        'pyscf_front',
        'pyscf_front/gui',
        'pyscf_front/core',
        'pyscf_front/database',
        'pyscf_front/utils',
        'pyscf_front/resources/styles'
    ]
    
    project_root = Path(__file__).parent
    results = []
    
    for dir_path in required_dirs:
        full_path = project_root / dir_path
        if full_path.exists():
            print(f"✓ Directory exists: {dir_path}")
            results.append(True)
        else:
            print(f"✗ Directory missing: {dir_path}")
            results.append(False)
    
    return all(results)

def test_gui_creation():
    """GUI作成テスト"""
    try:
        # QApplicationを作成（GUIなしモード）
        os.environ['QT_QPA_PLATFORM'] = 'offscreen'
        
        from PySide6.QtWidgets import QApplication, QLabel
        from PySide6.QtCore import Qt
        
        app = QApplication([])
        
        # 基本ウィジェット作成テスト
        label = QLabel("Test Label")
        label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        
        print("✓ Basic Qt widgets can be created")
        
        # アプリケーション設定テスト
        app.setApplicationName("Test App")
        app.setApplicationVersion("1.0.0")
        
        print("✓ Qt application settings work")
        
        return True
        
    except Exception as e:
        print(f"✗ GUI creation test failed: {e}")
        return False

def test_application_import():
    """アプリケーションインポートテスト"""
    try:
        # プロジェクトルートをパスに追加
        project_root = Path(__file__).parent
        if str(project_root) not in sys.path:
            sys.path.insert(0, str(project_root))
        
        # 設定クラスのインポートテスト
        from pyscf_front.utils.config import Config
        config = Config()
        print("✓ Config class imported and instantiated")
        
        # ロガーのインポートテスト
        from pyscf_front.utils.logger import setup_logger
        setup_logger()
        print("✓ Logger setup successful")
        
        # メインウィンドウのインポートテスト
        from pyscf_front.gui.main_window import MainWindow
        print("✓ MainWindow class imported successfully")
        
        return True
        
    except Exception as e:
        print(f"✗ Application import test failed: {e}")
        import traceback
        traceback.print_exc()
        return False

def run_all_tests():
    """全テストを実行"""
    print("=== PySCF_Front Development Environment Test ===\n")
    
    tests = [
        ("Python Version", test_python_version),
        ("PySide6 Import", test_pyside6),
        ("Scientific Packages", test_scientific_packages),
        ("Project Structure", test_project_structure),
        ("GUI Creation", test_gui_creation),
        ("Application Import", test_application_import),
    ]
    
    results = []
    for test_name, test_func in tests:
        print(f"\n--- {test_name} Test ---")
        try:
            result = test_func()
            results.append((test_name, result))
        except Exception as e:
            print(f"✗ {test_name} test crashed: {e}")
            results.append((test_name, False))
    
    print("\n=== Test Results Summary ===")
    passed = 0
    total = len(results)
    
    for test_name, result in results:
        status = "PASS" if result else "FAIL"
        print(f"{test_name}: {status}")
        if result:
            passed += 1
    
    print(f"\nTotal: {passed}/{total} tests passed")
    
    if passed == total:
        print("🎉 All tests passed! Development environment is ready.")
        return True
    else:
        print("❌ Some tests failed. Please check the setup.")
        return False

if __name__ == "__main__":
    success = run_all_tests()
    sys.exit(0 if success else 1)