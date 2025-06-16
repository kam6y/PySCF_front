#!/usr/bin/env python
"""開発用起動スクリプト"""
import sys
from pathlib import Path

# プロジェクトルートをパスに追加
project_root = Path(__file__).parent
sys.path.insert(0, str(project_root))

from PySide6.QtWidgets import QApplication
from PySide6.QtCore import Qt
from PySide6.QtGui import QIcon

def setup_application():
    """アプリケーションの基本設定"""
    app = QApplication(sys.argv)
    
    # 開発モード設定
    app.setOrganizationName("PySCF_Front")
    app.setOrganizationDomain("pyscf-front.org")
    app.setApplicationName("PySCF_Front Dev")
    app.setApplicationVersion("0.0.1-dev")
    
    # アイコンパスの設定
    icon_path = project_root / "pyscf_front" / "resources" / "icons"
    if icon_path.exists():
        app.setWindowIcon(QIcon(str(icon_path / "app_icon.png")))
    
    # スタイルシートの適用
    style_path = project_root / "pyscf_front" / "resources" / "styles" / "dark_theme.qss"
    if style_path.exists():
        with open(style_path, "r", encoding="utf-8") as f:
            app.setStyleSheet(f.read())
    
    # デバッグ情報
    print(f"Project root: {project_root}")
    print(f"Python version: {sys.version}")
    print(f"Qt version: {app.property('qt_version') or 'Unknown'}")
    
    return app

def main():
    """メイン関数"""
    try:
        app = setup_application()
        
        # メインウィンドウのインポートと作成
        try:
            from pyscf_front.gui.main_window import MainWindow
            window = MainWindow()
            window.show()
        except ImportError as e:
            print(f"Warning: Could not import MainWindow: {e}")
            print("Creating basic test window...")
            
            # 基本的なテストウィンドウを作成
            from PySide6.QtWidgets import QMainWindow, QLabel, QVBoxLayout, QWidget
            
            class TestWindow(QMainWindow):
                def __init__(self):
                    super().__init__()
                    self.setWindowTitle("PySCF_Front - Development Test")
                    self.setGeometry(100, 100, 800, 600)
                    
                    central_widget = QWidget()
                    layout = QVBoxLayout()
                    
                    label = QLabel("PySCF_Front Development Environment")
                    label.setAlignment(Qt.AlignmentFlag.AlignCenter)
                    layout.addWidget(label)
                    
                    central_widget.setLayout(layout)
                    self.setCentralWidget(central_widget)
            
            window = TestWindow()
            window.show()
        
        return app.exec()
        
    except Exception as e:
        print(f"Error starting application: {e}")
        import traceback
        traceback.print_exc()
        return 1

if __name__ == "__main__":
    exit_code = main()
    sys.exit(exit_code)