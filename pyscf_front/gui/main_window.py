"""
PySCF_Front メインウィンドウ
"""
from PySide6.QtWidgets import (
    QMainWindow, QWidget, QVBoxLayout, QDockWidget,
    QLabel, QTextEdit, QPushButton
)
from PySide6.QtCore import Qt, Signal, Slot
from PySide6.QtGui import QAction, QKeySequence
from loguru import logger

from pyscf_front.gui.widgets.molecule_viewer import MoleculeManagerWidget
from pyscf_front.core.molecule import MoleculeManager
from pyscf_front.core.calculation_engine_unified import UnifiedCalculationEngine
from pyscf_front.gui.dialogs.calculation_dialog import CalculationSettingsDialog

class MainWindow(QMainWindow):
    """PySCF_Front メインウィンドウ"""
    
    # シグナル定義
    calculation_started = Signal(str)  # calculation_id
    calculation_completed = Signal(str, dict)  # calculation_id, results
    
    def __init__(self):
        super().__init__()
        logger.info("Initializing MainWindow")
        
        # 分子管理システム
        self.molecule_manager = MoleculeManager()
        
        # 計算エンジン（統合版・SQLiteデータベース専用）
        try:
            self.calculation_engine = UnifiedCalculationEngine()
        except RuntimeError as e:
            logger.error(f"Failed to initialize calculation engine: {e}")
            self.show_database_error(str(e))
            return
        
        # ウィンドウ基本設定
        self.setWindowTitle("PySCF_Front v0.0.1")
        self.setGeometry(100, 100, 1400, 900)
        
        # 中央ウィジェット設定
        self.setup_central_widget()
        
        # メニューバーとツールバー
        self.setup_menu_bar()
        self.setup_status_bar()
        
        # ドックウィジェット
        self.setup_dock_widgets()
        
        # シグナル接続
        self.connect_signals()
        self.connect_calculation_signals()
        
        # 中断された計算を再開
        self.resume_interrupted_calculations()
        
        logger.info("MainWindow initialization completed")
    
    def setup_central_widget(self):
        """中央ウィジェットの設定"""
        # 分子管理ウィジェットを中央に配置
        self.molecule_widget = MoleculeManagerWidget()
        self.setCentralWidget(self.molecule_widget)
    
    def setup_menu_bar(self):
        """メニューバーの設定"""
        menubar = self.menuBar()
        
        # ファイルメニュー
        file_menu = menubar.addMenu("ファイル(&F)")
        
        new_action = QAction("新規プロジェクト(&N)", self)
        new_action.setShortcut(QKeySequence.StandardKey.New)
        new_action.setStatusTip("新しいプロジェクトを作成")
        new_action.triggered.connect(self.new_project)
        file_menu.addAction(new_action)
        
        open_action = QAction("プロジェクトを開く(&O)", self)
        open_action.setShortcut(QKeySequence.StandardKey.Open)
        open_action.setStatusTip("既存のプロジェクトを開く")
        open_action.triggered.connect(self.open_project)
        file_menu.addAction(open_action)
        
        file_menu.addSeparator()
        
        exit_action = QAction("終了(&X)", self)
        exit_action.setShortcut(QKeySequence.StandardKey.Quit)
        exit_action.setStatusTip("アプリケーションを終了")
        exit_action.triggered.connect(self.close)
        file_menu.addAction(exit_action)
        
        # 計算メニュー
        calc_menu = menubar.addMenu("計算(&C)")
        
        start_calc_action = QAction("計算設定(&S)", self)
        start_calc_action.setShortcut("Ctrl+R")
        start_calc_action.setStatusTip("計算設定ダイアログを開く")
        start_calc_action.triggered.connect(self.show_calculation_dialog)
        calc_menu.addAction(start_calc_action)
        
        quick_calc_action = QAction("クイック計算(&Q)", self)
        quick_calc_action.setStatusTip("デフォルト設定で計算を開始")
        quick_calc_action.triggered.connect(self.start_quick_calculation)
        calc_menu.addAction(quick_calc_action)
        
        calc_menu.addSeparator()
        
        job_list_action = QAction("ジョブリスト(&L)", self)
        job_list_action.setStatusTip("計算ジョブの一覧を表示")
        job_list_action.triggered.connect(self.show_job_list)
        calc_menu.addAction(job_list_action)
        
        calc_menu.addSeparator()
        
        history_action = QAction("計算履歴(&H)", self)
        history_action.setStatusTip("過去の計算履歴を表示")
        history_action.triggered.connect(self.show_calculation_history)
        calc_menu.addAction(history_action)
        
        instances_action = QAction("インスタンス管理(&I)", self)
        instances_action.setStatusTip("プロジェクトインスタンスを管理")
        instances_action.triggered.connect(self.show_instances_manager)
        calc_menu.addAction(instances_action)
        
        # ツールメニュー
        tools_menu = menubar.addMenu("ツール(&T)")
        
        settings_action = QAction("設定(&S)", self)
        settings_action.setStatusTip("アプリケーション設定")
        settings_action.triggered.connect(self.show_settings)
        tools_menu.addAction(settings_action)
        
        # ヘルプメニュー
        help_menu = menubar.addMenu("ヘルプ(&H)")
        
        about_action = QAction("PySCF_Frontについて(&A)", self)
        about_action.setStatusTip("このアプリケーションについて")
        about_action.triggered.connect(self.show_about)
        help_menu.addAction(about_action)
    
    def setup_status_bar(self):
        """ステータスバーの設定"""
        self.status_bar = self.statusBar()
        self.status_bar.showMessage("準備完了")
        
        # ステータス情報ラベル
        self.status_label = QLabel("PySCF_Front v2.0")
        self.status_bar.addPermanentWidget(self.status_label)
    
    def setup_dock_widgets(self):
        """ドックウィジェットの設定"""
        # 計算設定パネル
        self.calc_dock = QDockWidget("計算設定", self)
        self.calc_widget = QWidget()
        calc_layout = QVBoxLayout()
        
        calc_layout.addWidget(QLabel("計算手法設定"))
        
        # 計算設定テキスト
        self.calc_settings_text = QTextEdit()
        self.calc_settings_text.setPlainText("手法: B3LYP\n基底関数: 6-31G(d)\n電荷: 0\nスピン多重度: 1")
        self.calc_settings_text.setMaximumHeight(150)
        calc_layout.addWidget(self.calc_settings_text)
        
        # 計算ボタン
        from PySide6.QtWidgets import QHBoxLayout
        calc_btn_layout = QHBoxLayout()
        self.calc_settings_btn = QPushButton("計算設定")
        self.quick_calc_btn = QPushButton("クイック計算")
        
        self.calc_settings_btn.clicked.connect(self.show_calculation_dialog)
        self.quick_calc_btn.clicked.connect(self.start_quick_calculation)
        
        calc_btn_layout.addWidget(self.calc_settings_btn)
        calc_btn_layout.addWidget(self.quick_calc_btn)
        calc_layout.addLayout(calc_btn_layout)
        
        calc_layout.addStretch()
        self.calc_widget.setLayout(calc_layout)
        self.calc_dock.setWidget(self.calc_widget)
        self.addDockWidget(Qt.DockWidgetArea.RightDockWidgetArea, self.calc_dock)
        
        # ジョブモニター
        self.job_dock = QDockWidget("ジョブモニター", self)
        self.job_widget = QTextEdit()
        self.job_widget.setPlainText("実行中のジョブ\n（現在なし）")
        self.job_widget.setMaximumHeight(200)
        self.job_dock.setWidget(self.job_widget)
        self.addDockWidget(Qt.DockWidgetArea.BottomDockWidgetArea, self.job_dock)
    
    def connect_signals(self):
        """シグナル接続"""
        # 分子選択シグナル
        self.molecule_widget.molecule_selected.connect(self.on_molecule_selected)
        
        # 計算関連シグナル
        self.calculation_started.connect(self.on_calculation_started)
        self.calculation_completed.connect(self.on_calculation_completed)
    
    def connect_calculation_signals(self):
        """計算エンジンのシグナルを接続"""
        # 計算エンジンからの進捗シグナルを接続
        if hasattr(self.calculation_engine, 'signals'):
            self.calculation_engine.signals.started.connect(self.on_calculation_started)
            self.calculation_engine.signals.progress.connect(self.on_calculation_progress)
            self.calculation_engine.signals.completed.connect(self.on_calculation_completed)
            self.calculation_engine.signals.failed.connect(self.on_calculation_failed)
    
    def resume_interrupted_calculations(self):
        """中断された計算を再開"""
        try:
            resumed_jobs = self.calculation_engine.resume_interrupted_calculations()
            if resumed_jobs:
                self.status_bar.showMessage(f"{len(resumed_jobs)}件の中断された計算を再開しました", 5000)
                logger.info(f"Resumed {len(resumed_jobs)} interrupted calculations")
                self.update_job_monitor()
        except Exception as e:
            logger.error(f"Failed to resume interrupted calculations: {e}")
            self.status_bar.showMessage("中断された計算の再開に失敗しました", 3000)
    
    @Slot(object)
    def on_molecule_selected(self, molecule):
        """分子選択時の処理"""
        logger.info(f"Molecule selected: {molecule.name}")
        self.status_bar.showMessage(f"分子を選択: {molecule.name}")
        
        # 計算設定を更新
        settings_text = f"手法: B3LYP\n基底関数: 6-31G(d)\n"
        settings_text += f"電荷: {molecule.charge}\n"
        settings_text += f"スピン多重度: {molecule.multiplicity}\n"
        settings_text += f"\n分子: {molecule.name}\n"
        settings_text += f"分子式: {molecule.get_molecular_formula()}\n"
        settings_text += f"原子数: {len(molecule.atoms)}"
        
        self.calc_settings_text.setPlainText(settings_text)
    
    
    @Slot()
    def new_project(self):
        """新規プロジェクト作成"""
        logger.info("New project requested")
        self.status_bar.showMessage("新規プロジェクト作成（未実装）", 2000)
    
    @Slot()
    def open_project(self):
        """プロジェクトを開く"""
        logger.info("Open project requested")
        self.status_bar.showMessage("プロジェクトを開く（未実装）", 2000)
    
    @Slot()
    def show_calculation_dialog(self):
        """計算設定ダイアログを表示"""
        if not hasattr(self.molecule_widget, 'current_molecule') or not self.molecule_widget.current_molecule:
            self.status_bar.showMessage("分子を選択してください", 3000)
            return
        
        molecule = self.molecule_widget.current_molecule
        dialog = CalculationSettingsDialog(molecule, self.calculation_engine, self)
        dialog.exec()
    
    @Slot()
    def start_quick_calculation(self):
        """クイック計算開始"""
        if not hasattr(self.molecule_widget, 'current_molecule') or not self.molecule_widget.current_molecule:
            self.status_bar.showMessage("分子を選択してください", 3000)
            return
        
        molecule = self.molecule_widget.current_molecule
        
        try:
            # SQLiteデータベースへの永続化計算
            job_id = self.calculation_engine.submit_calculation(
                molecule,
                method="B3LYP",
                basis_set="6-31G(d)",
                instance_name=f"QuickCalc_{molecule.name}"
            )
            self.status_bar.showMessage(f"クイック計算を開始: {molecule.name} (Job: {job_id[:8]})", 5000)
            logger.info(f"Started quick calculation for {molecule.name}: {job_id}")
            
        except Exception as e:
            self.status_bar.showMessage(f"計算開始エラー: {str(e)}", 5000)
            logger.error(f"Failed to start quick calculation: {e}")
    
    @Slot()
    def show_job_list(self):
        """ジョブリストを表示"""
        jobs = self.calculation_engine.list_jobs()
        
        if not jobs:
            self.status_bar.showMessage("実行中のジョブはありません", 3000)
            return
        
        # 簡単なジョブリスト表示
        job_text = "計算ジョブ一覧:\n\n"
        for job in jobs:
            job_text += f"ID: {job['id'][:8]}\n"
            job_text += f"分子: {job['molecule_name']}\n"
            job_text += f"手法: {job['method']}/{job['basis_set']}\n"
            job_text += f"状態: {job['status']}\n"
            if job['progress'] > 0:
                job_text += f"進捗: {job['progress']:.1%}\n"
            job_text += "\n"
        
        self.job_widget.setPlainText(job_text)
    
    @Slot()
    def show_settings(self):
        """設定ダイアログ表示"""
        logger.info("Settings dialog requested")
        self.status_bar.showMessage("設定ダイアログ（未実装）", 2000)
    
    @Slot()
    def show_about(self):
        """アバウトダイアログ表示"""
        from PySide6.QtWidgets import QMessageBox
        
        QMessageBox.about(
            self,
            "PySCF_Frontについて",
            """
            <h3>PySCF_Front v0.0.1</h3>
            <p>PySCFライブラリを使用した量子化学計算のためのGUIアプリケーション</p>
            <p>開発環境: PySide6 + Python</p>
            <p>実装済み機能:</p>
            <ul>
            <li>分子構造入力 (SMILES, XYZ, 手動, プリセット)</li>
            <li>分子表示・管理</li>
            <li>簡易3D可視化</li>
            <li>PySCF量子化学計算エンジン</li>
            <li>計算進捗監視とジョブ管理</li>
            </ul>
            <p>ライセンス: MIT License</p>
            """
        )
    
    @Slot(str)
    def on_calculation_started(self, job_id):
        """計算開始時の処理"""
        logger.info(f"Calculation started: {job_id}")
        self.status_bar.showMessage(f"計算開始: {job_id[:8]}")
        
        # ジョブモニターを更新
        job = self.calculation_engine.get_job_status(job_id)
        if job:
            self.update_job_monitor()
    
    @Slot(str, float, str)
    def on_calculation_progress(self, job_id, progress, message):
        """計算進捗更新時の処理"""
        logger.info(f"Calculation progress {job_id[:8]}: {progress:.1%} - {message}")
        self.status_bar.showMessage(f"計算進捗 {job_id[:8]}: {progress:.1%} - {message}")
        
        # ジョブモニターを更新
        self.update_job_monitor()
    
    @Slot(str, dict)
    def on_calculation_completed(self, job_id, results):
        """計算完了時の処理"""
        logger.info(f"Calculation completed: {job_id}")
        
        job = self.calculation_engine.get_job_status(job_id)
        if job:
            energy = results.get('total_energy', 0.0)
            converged = results.get('converged', False)
            
            status_msg = f"計算完了: {job.molecule.name} E={energy:.6f} Ha"
            if not converged:
                status_msg += " (収束せず)"
            
            self.status_bar.showMessage(status_msg, 10000)
            
            # 結果をジョブモニターに表示
            self.update_job_monitor()
            self.show_calculation_results(job_id, results)
    
    @Slot(str, str)
    def on_calculation_failed(self, job_id, error_msg):
        """計算失敗時の処理"""
        logger.error(f"Calculation failed {job_id}: {error_msg}")
        self.status_bar.showMessage(f"計算失敗 {job_id[:8]}: {error_msg}", 10000)
        
        # ジョブモニターを更新
        self.update_job_monitor()
    
    def update_job_monitor(self):
        """ジョブモニターを更新"""
        jobs = self.calculation_engine.list_jobs()
        
        if not jobs:
            self.job_widget.setPlainText("実行中のジョブ\n（現在なし）")
            return
        
        job_text = "計算ジョブ状況:\n\n"
        
        for job in jobs[-5:]:  # 最新5件のみ表示
            status_icon = {
                'pending': '⏳',
                'running': '🔄',
                'completed': '✅',
                'failed': '❌',
                'cancelled': '🚫'
            }.get(job['status'], '❓')
            
            job_text += f"{status_icon} {job['id'][:8]} - {job['molecule_name']}\n"
            job_text += f"   {job['method']}/{job['basis_set']}\n"
            
            if job['status'] == 'running':
                job_text += f"   進捗: {job['progress']:.1%}\n"
            elif job['status'] == 'completed':
                job_text += f"   完了\n"
            elif job['status'] == 'failed':
                job_text += f"   失敗\n"
            
            job_text += "\n"
        
        self.job_widget.setPlainText(job_text)
    
    @Slot()
    def show_calculation_history(self):
        """計算履歴を表示"""
        try:
            history = self.calculation_engine.get_calculation_history()
            
            if not history:
                self.status_bar.showMessage("計算履歴がありません", 3000)
                return
            
            from PySide6.QtWidgets import QDialog, QVBoxLayout, QTextEdit, QPushButton
            from PySide6.QtGui import QFont
            
            dialog = QDialog(self)
            dialog.setWindowTitle("計算履歴")
            dialog.setMinimumSize(800, 600)
            
            layout = QVBoxLayout()
            
            text_edit = QTextEdit()
            text_edit.setReadOnly(True)
            text_edit.setFont(QFont("Courier", 9))
            
            history_text = "計算履歴:\n\n"
            for i, calc in enumerate(history, 1):
                history_text += f"=== 計算 {i} ===\n"
                history_text += f"インスタンス: {calc.get('instance_name', 'Unknown')}\n"
                history_text += f"分子: {calc.get('molecule_name', 'Unknown')}\n"
                history_text += f"手法: {calc.get('method', 'Unknown')}/{calc.get('basis_set', 'Unknown')}\n"
                history_text += f"状態: {calc.get('status', 'Unknown')}\n"
                history_text += f"作成日時: {calc.get('created_at', 'Unknown')}\n"
                if calc.get('completed_at'):
                    history_text += f"完了日時: {calc['completed_at']}\n"
                if calc.get('total_energy') is not None:
                    try:
                        energy_val = float(calc['total_energy'])
                        history_text += f"エネルギー: {energy_val:.8f} Ha\n"
                    except (ValueError, TypeError):
                        history_text += f"エネルギー: {calc['total_energy']}\n"
                if calc.get('error_message'):
                    history_text += f"エラー: {calc['error_message']}\n"
                history_text += "-" * 50 + "\n\n"
            
            text_edit.setPlainText(history_text)
            layout.addWidget(text_edit)
            
            close_btn = QPushButton("閉じる")
            close_btn.clicked.connect(dialog.close)
            layout.addWidget(close_btn)
            
            dialog.setLayout(layout)
            dialog.exec()
            
        except Exception as e:
            self.status_bar.showMessage(f"履歴表示エラー: {str(e)}", 5000)
            logger.error(f"Failed to show calculation history: {e}")
    
    @Slot()
    def show_instances_manager(self):
        """インスタンス管理ダイアログを表示"""
        try:
            instances = self.calculation_engine.list_all_instances()
            
            from PySide6.QtWidgets import (QDialog, QVBoxLayout, QHBoxLayout, 
                                         QListWidget, QListWidgetItem, QPushButton, 
                                         QTextEdit, QMessageBox, QInputDialog)
            
            dialog = QDialog(self)
            dialog.setWindowTitle("インスタンス管理")
            dialog.setMinimumSize(900, 700)
            
            layout = QHBoxLayout()
            
            # 左側: インスタンスリスト
            left_layout = QVBoxLayout()
            left_layout.addWidget(QLabel("プロジェクトインスタンス"))
            
            instance_list = QListWidget()
            for instance in instances:
                item_text = f"{instance['name']} ({instance.get('status', 'unknown')})"
                item = QListWidgetItem(item_text)
                item.setData(Qt.ItemDataRole.UserRole, instance['id'])
                instance_list.addItem(item)
            
            left_layout.addWidget(instance_list)
            
            # ボタン
            btn_layout = QHBoxLayout()
            view_btn = QPushButton("詳細表示")
            delete_btn = QPushButton("削除")
            refresh_btn = QPushButton("更新")
            
            btn_layout.addWidget(view_btn)
            btn_layout.addWidget(delete_btn)
            btn_layout.addWidget(refresh_btn)
            left_layout.addLayout(btn_layout)
            
            left_widget = QWidget()
            left_widget.setLayout(left_layout)
            left_widget.setMaximumWidth(300)
            
            # 右側: 詳細情報
            right_layout = QVBoxLayout()
            right_layout.addWidget(QLabel("インスタンス詳細"))
            
            details_text = QTextEdit()
            details_text.setReadOnly(True)
            details_text.setFont(QFont("Courier", 9))
            details_text.setPlainText("インスタンスを選択してください...")
            
            right_layout.addWidget(details_text)
            
            right_widget = QWidget()
            right_widget.setLayout(right_layout)
            
            layout.addWidget(left_widget)
            layout.addWidget(right_widget)
            
            # シグナル接続
            def on_instance_selected():
                current_item = instance_list.currentItem()
                if current_item:
                    instance_id = current_item.data(Qt.ItemDataRole.UserRole)
                    try:
                        instance_details = self.calculation_engine.get_instance_details(instance_id)
                        if instance_details:
                            details_text.setPlainText(f"""インスタンス詳細:

ID: {instance_details['id']}
名前: {instance_details['name']}
説明: {instance_details.get('description', 'N/A')}
作成日時: {instance_details.get('created_at', 'N/A')}
状態: {instance_details.get('status', 'N/A')}

分子情報:
{instance_details.get('molecule_info', '分子情報なし')}

計算一覧:
{instance_details.get('calculations_summary', '計算なし')}

結果サマリー:
{instance_details.get('results_summary', '結果なし')}
""")
                    except Exception as e:
                        details_text.setPlainText(f"詳細取得エラー: {str(e)}")
            
            def delete_instance():
                current_item = instance_list.currentItem()
                if current_item:
                    instance_id = current_item.data(Qt.ItemDataRole.UserRole)
                    instance_name = current_item.text().split(' (')[0]
                    
                    reply = QMessageBox.question(
                        dialog,
                        "削除確認",
                        f"インスタンス '{instance_name}' を削除しますか？\n\n"
                        "この操作は取り消せません。",
                        QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No
                    )
                    
                    if reply == QMessageBox.StandardButton.Yes:
                        try:
                            if self.calculation_engine.delete_instance(instance_id):
                                QMessageBox.information(dialog, "削除完了", "インスタンスを削除しました")
                                # リストを更新
                                instance_list.takeItem(instance_list.currentRow())
                                details_text.clear()
                            else:
                                QMessageBox.warning(dialog, "削除失敗", "インスタンスの削除に失敗しました")
                        except Exception as e:
                            QMessageBox.critical(dialog, "エラー", f"削除中にエラーが発生しました:\n{str(e)}")
            
            def refresh_list():
                instance_list.clear()
                try:
                    updated_instances = self.calculation_engine.list_all_instances()
                    for instance in updated_instances:
                        item_text = f"{instance['name']} ({instance.get('status', 'unknown')})"
                        item = QListWidgetItem(item_text)
                        item.setData(Qt.ItemDataRole.UserRole, instance['id'])
                        instance_list.addItem(item)
                except Exception as e:
                    QMessageBox.warning(dialog, "更新失敗", f"リストの更新に失敗しました:\n{str(e)}")
            
            instance_list.itemSelectionChanged.connect(on_instance_selected)
            view_btn.clicked.connect(on_instance_selected)
            delete_btn.clicked.connect(delete_instance)
            refresh_btn.clicked.connect(refresh_list)
            
            dialog.setLayout(layout)
            dialog.exec()
            
        except Exception as e:
            self.status_bar.showMessage(f"インスタンス管理エラー: {str(e)}", 5000)
            logger.error(f"Failed to show instances manager: {e}")
    
    def show_calculation_results(self, job_id: str, results: dict):
        """計算結果を表示"""
        from PySide6.QtWidgets import QMessageBox
        
        job = self.calculation_engine.get_job_status(job_id)
        if not job:
            return
        
        # 結果表示用テキスト作成
        result_text = f"""計算結果 - {job.molecule.name}

基本情報:
  計算手法: {job.method}
  基底関数: {job.basis_set}
  計算時間: {results.get('calculation_time', 0):.2f}秒
  収束: {'はい' if results.get('converged', False) else 'いいえ'}

エネルギー:
  全エネルギー: {results.get('total_energy', 0):.8f} Hartree
  HOMO エネルギー: {results.get('homo_energy', 0):.6f} Hartree
  LUMO エネルギー: {results.get('lumo_energy', 0):.6f} Hartree
  HOMO-LUMO ギャップ: {results.get('homo_lumo_gap', 0):.6f} Hartree

分子特性:
  双極子モーメント: {results.get('dipole_moment', 0):.4f} a.u.
  """
        
        if results.get('atomic_charges'):
            result_text += "\nMulliken電荷:\n"
            for i, charge in enumerate(results['atomic_charges']):
                atom_symbol = job.molecule.atoms[i].symbol
                result_text += f"  {atom_symbol}{i+1}: {charge:+.4f}\n"
        
        # メッセージボックスで表示
        msg = QMessageBox(self)
        msg.setWindowTitle(f"計算結果 - {job.molecule.name}")
        msg.setText(result_text)
        msg.setIcon(QMessageBox.Icon.Information)
        msg.exec()
    
    def show_database_error(self, error_message: str):
        """データベースエラーの表示"""
        from PySide6.QtWidgets import QMessageBox
        
        QMessageBox.critical(
            self,
            "データベースエラー",
            f"SQLiteデータベースの初期化に失敗しました:\n\n{error_message}\n\n"
            "アプリケーションを終了します。"
        )
        self.close()
    
    def closeEvent(self, event):
        """アプリケーション終了時の処理"""
        logger.info("Application closing")
        event.accept()