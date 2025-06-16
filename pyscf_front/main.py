"""
PySCF_Front メインアプリケーション
"""
import sys
from pathlib import Path
from PySide6.QtWidgets import QApplication
from PySide6.QtCore import QTranslator, QLocale, Qt
from PySide6.QtGui import QIcon
from loguru import logger

# プロジェクトルートをパスに追加
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

from pyscf_front.utils.config import Config
from pyscf_front.utils.logger import setup_logger


class PySCFFrontApplication(QApplication):
    """PySCF_Front アプリケーションクラス"""
    
    def __init__(self, argv):
        super().__init__(argv)
        self.config = Config()
        self.setup_application()
        self.setup_translation()
        self.setup_style()
        
    def setup_application(self):
        """アプリケーション基本設定"""
        self.setOrganizationName("PySCF_Front")
        self.setOrganizationDomain("pyscf-front.org")
        self.setApplicationName("PySCF_Front")
        self.setApplicationVersion("0.0.1")
        
        # アイコン設定
        icon_path = project_root / "pyscf_front" / "resources" / "icons" / "app_icon.png"
        if icon_path.exists():
            self.setWindowIcon(QIcon(str(icon_path)))
            
    def setup_translation(self):
        """国際化設定"""
        translator = QTranslator()
        locale = self.config.get('language', QLocale.system().name())
        
        translation_path = project_root / "pyscf_front" / "translations" / f"pyscf_front_{locale}.qm"
        if translation_path.exists() and translator.load(str(translation_path)):
            self.installTranslator(translator)
            logger.info(f"Loaded translation: {locale}")
            
    def setup_style(self):
        """スタイル設定"""
        theme = self.config.get('theme', 'dark')
        style_path = project_root / "pyscf_front" / "resources" / "styles" / f"{theme}_theme.qss"
        
        if style_path.exists():
            with open(style_path, "r", encoding="utf-8") as f:
                self.setStyleSheet(f.read())
            logger.info(f"Applied theme: {theme}")


def main():
    """メイン関数"""
    try:
        # ロガー設定
        setup_logger()
        logger.info("Starting PySCF_Front...")
        
        # Qt アプリケーション作成
        app = PySCFFrontApplication(sys.argv)
        
        # メインウィンドウ作成
        try:
            from pyscf_front.gui.main_window import MainWindow
            window = MainWindow()
            window.show()
            logger.info("Main window created and shown")
        except ImportError as e:
            logger.warning(f"Could not import MainWindow: {e}")
            logger.info("Running in minimal mode...")
            
            # 最小限のテストウィンドウ
            from PySide6.QtWidgets import QMainWindow, QLabel, QVBoxLayout, QWidget
            
            class MinimalWindow(QMainWindow):
                def __init__(self):
                    super().__init__()
                    self.setWindowTitle("PySCF_Front - Minimal Mode")
                    self.setGeometry(200, 200, 600, 400)
                    
                    central_widget = QWidget()
                    layout = QVBoxLayout()
                    
                    label = QLabel("PySCF_Front - Environment Ready")
                    label.setAlignment(Qt.AlignmentFlag.AlignCenter)
                    layout.addWidget(label)
                    
                    status_label = QLabel("Main application components not yet implemented")
                    status_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
                    layout.addWidget(status_label)
                    
                    central_widget.setLayout(layout)
                    self.setCentralWidget(central_widget)
            
            window = MinimalWindow()
            window.show()
        
        logger.info("Application started successfully")
        return app.exec()
        
    except Exception as e:
        logger.error(f"Error starting application: {e}")
        logger.exception("Full traceback:")
        return 1


if __name__ == "__main__":
    exit_code = main()
    sys.exit(exit_code)