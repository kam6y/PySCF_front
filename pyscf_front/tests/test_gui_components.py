"""
GUIコンポーネントのテストスイート
"""
import pytest
from unittest.mock import Mock, patch, MagicMock
from PySide6.QtWidgets import QApplication, QWidget
from PySide6.QtCore import QTimer, Qt
from PySide6.QtTest import QTest

# Qt テスト用の設定
@pytest.fixture(scope="session")
def qapp():
    """QApplication のセッションフィクスチャ"""
    app = QApplication.instance()
    if app is None:
        app = QApplication([])
    yield app
    # app.quit() # テスト終了時にアプリケーションを終了


class TestMainWindow:
    """MainWindow の GUI テスト"""
    
    def test_main_window_initialization(self, qapp, mock_pyscf):
        """メインウィンドウの初期化テスト"""
        with patch('pyscf_front.gui.main_window.UnifiedCalculationEngine') as mock_engine_cls, \
             patch('pyscf_front.gui.main_window.MoleculeManager') as mock_mol_mgr_cls:
            
            # モックの設定
            mock_engine = Mock()
            mock_engine.resume_interrupted_calculations.return_value = []
            mock_engine.signals = Mock()
            mock_engine.signals.started = Mock()
            mock_engine.signals.progress = Mock()
            mock_engine.signals.completed = Mock()
            mock_engine.signals.failed = Mock()
            
            mock_engine_cls.return_value = mock_engine
            mock_mol_mgr_cls.return_value = Mock()
            
            from pyscf_front.gui.main_window import MainWindow
            
            window = MainWindow()
            
            # 基本的な属性の確認
            assert window.windowTitle() == "PySCF_Front v0.0.1"
            assert hasattr(window, 'molecule_manager')
            assert hasattr(window, 'calculation_engine')
            assert hasattr(window, 'molecule_widget')
            
            # モックが呼ばれたことを確認
            mock_engine_cls.assert_called_once()
            mock_engine.resume_interrupted_calculations.assert_called_once()
            
            window.close()
    
    def test_menu_bar_creation(self, qapp, mock_pyscf):
        """メニューバーの作成テスト"""
        with patch('pyscf_front.gui.main_window.UnifiedCalculationEngine') as mock_engine_cls, \
             patch('pyscf_front.gui.main_window.MoleculeManager') as mock_mol_mgr_cls:
            
            mock_engine = Mock()
            mock_engine.resume_interrupted_calculations.return_value = []
            mock_engine.signals = Mock()
            mock_engine.signals.started = Mock()
            mock_engine.signals.progress = Mock()
            mock_engine.signals.completed = Mock()
            mock_engine.signals.failed = Mock()
            
            mock_engine_cls.return_value = mock_engine
            mock_mol_mgr_cls.return_value = Mock()
            
            from pyscf_front.gui.main_window import MainWindow
            
            window = MainWindow()
            
            # メニューバーの存在確認
            menu_bar = window.menuBar()
            assert menu_bar is not None
            
            # メニューの確認
            menus = [action.text() for action in menu_bar.actions()]
            assert "ファイル(&F)" in menus
            assert "計算(&C)" in menus
            assert "ツール(&T)" in menus
            assert "ヘルプ(&H)" in menus
            
            window.close()
    
    def test_dock_widgets_creation(self, qapp, mock_pyscf):
        """ドックウィジェットの作成テスト"""
        with patch('pyscf_front.gui.main_window.UnifiedCalculationEngine') as mock_engine_cls, \
             patch('pyscf_front.gui.main_window.MoleculeManager') as mock_mol_mgr_cls:
            
            mock_engine = Mock()
            mock_engine.resume_interrupted_calculations.return_value = []
            mock_engine.signals = Mock()
            mock_engine.signals.started = Mock()
            mock_engine.signals.progress = Mock()
            mock_engine.signals.completed = Mock()
            mock_engine.signals.failed = Mock()
            
            mock_engine_cls.return_value = mock_engine
            mock_mol_mgr_cls.return_value = Mock()
            
            from pyscf_front.gui.main_window import MainWindow
            
            window = MainWindow()
            
            # ドックウィジェットの確認
            assert hasattr(window, 'calc_dock')
            assert hasattr(window, 'job_dock')
            assert window.calc_dock.windowTitle() == "計算設定"
            assert window.job_dock.windowTitle() == "ジョブモニター"
            
            window.close()
    
    def test_quick_calculation_button_click(self, qapp, qtbot, sample_molecule, mock_pyscf):
        """クイック計算ボタンのクリックテスト"""
        with patch('pyscf_front.gui.main_window.UnifiedCalculationEngine') as mock_engine_cls, \
             patch('pyscf_front.gui.main_window.MoleculeManager') as mock_mol_mgr_cls:
            
            mock_engine = Mock()
            mock_engine.resume_interrupted_calculations.return_value = []
            mock_engine.submit_calculation_with_persistence.return_value = "test-job-id"
            mock_engine.signals = Mock()
            mock_engine.signals.started = Mock()
            mock_engine.signals.progress = Mock()
            mock_engine.signals.completed = Mock()
            mock_engine.signals.failed = Mock()
            
            mock_engine_cls.return_value = mock_engine
            mock_mol_mgr_cls.return_value = Mock()
            
            from pyscf_front.gui.main_window import MainWindow
            
            window = MainWindow()
            qtbot.addWidget(window)
            
            # 分子を設定
            window.molecule_widget.current_molecule = sample_molecule
            
            # クイック計算ボタンをクリック
            qtbot.mouseClick(window.quick_calc_btn, Qt.MouseButton.LeftButton)
            
            # 計算が投入されたことを確認
            mock_engine.submit_calculation_with_persistence.assert_called_once()
            args, kwargs = mock_engine.submit_calculation_with_persistence.call_args
            assert kwargs['method'] == "B3LYP"
            assert kwargs['basis_set'] == "6-31G(d)"
            
            window.close()
    
    def test_molecule_selection_updates_display(self, qapp, qtbot, sample_molecule, mock_pyscf):
        """分子選択時の表示更新テスト"""
        with patch('pyscf_front.gui.main_window.UnifiedCalculationEngine') as mock_engine_cls, \
             patch('pyscf_front.gui.main_window.MoleculeManager') as mock_mol_mgr_cls:
            
            mock_engine = Mock()
            mock_engine.resume_interrupted_calculations.return_value = []
            mock_engine.signals = Mock()
            mock_engine.signals.started = Mock()
            mock_engine.signals.progress = Mock()
            mock_engine.signals.completed = Mock()
            mock_engine.signals.failed = Mock()
            
            mock_engine_cls.return_value = mock_engine
            mock_mol_mgr_cls.return_value = Mock()
            
            from pyscf_front.gui.main_window import MainWindow
            
            window = MainWindow()
            qtbot.addWidget(window)
            
            # 分子選択をシミュレート
            window.on_molecule_selected(sample_molecule)
            
            # 計算設定テキストが更新されたことを確認
            settings_text = window.calc_settings_text.toPlainText()
            assert sample_molecule.name in settings_text
            assert f"電荷: {sample_molecule.charge}" in settings_text
            assert f"スピン多重度: {sample_molecule.multiplicity}" in settings_text
            
            window.close()


class TestMoleculeViewer:
    """MoleculeViewer の GUI テスト"""
    
    def test_molecule_info_widget_initialization(self, qapp):
        """MoleculeInfoWidget の初期化テスト"""
        from pyscf_front.gui.widgets.molecule_viewer import MoleculeInfoWidget
        
        widget = MoleculeInfoWidget()
        
        # 基本的な属性の確認
        assert hasattr(widget, 'name_label')
        assert hasattr(widget, 'formula_label')
        assert hasattr(widget, 'atoms_label')
        assert hasattr(widget, 'charge_label')
        assert hasattr(widget, 'multiplicity_label')
        assert hasattr(widget, 'coords_text')
        
        # 初期状態の確認
        assert widget.name_label.text() == "名前: -"
        assert widget.formula_label.text() == "分子式: -"
        
        widget.close()
    
    def test_molecule_info_display(self, qapp, sample_molecule):
        """分子情報の表示テスト"""
        from pyscf_front.gui.widgets.molecule_viewer import MoleculeInfoWidget
        
        widget = MoleculeInfoWidget()
        
        # 分子を設定
        widget.set_molecule(sample_molecule)
        
        # 表示内容の確認
        assert widget.name_label.text() == f"名前: {sample_molecule.name}"
        assert widget.formula_label.text() == f"分子式: {sample_molecule.get_molecular_formula()}"
        assert widget.atoms_label.text() == f"原子数: {len(sample_molecule.atoms)}"
        assert widget.charge_label.text() == f"電荷: {sample_molecule.charge}"
        assert widget.multiplicity_label.text() == f"スピン多重度: {sample_molecule.multiplicity}"
        
        # 座標テキストの確認
        coords_text = widget.coords_text.toPlainText()
        assert "原子" in coords_text
        assert "X座標" in coords_text
        assert "Y座標" in coords_text
        assert "Z座標" in coords_text
        
        widget.close()
    
    def test_molecule_list_widget_functionality(self, qapp, qtbot, sample_molecule, sample_molecule_methane):
        """MoleculeListWidget の機能テスト"""
        from pyscf_front.gui.widgets.molecule_viewer import MoleculeListWidget
        
        widget = MoleculeListWidget()
        qtbot.addWidget(widget)
        
        # 分子を追加
        widget.add_molecule(sample_molecule)
        widget.add_molecule(sample_molecule_methane)
        
        # リストアイテムが追加されたことを確認
        assert widget.molecule_list.count() == 2
        
        # 分子IDが正しく保存されていることを確認
        assert sample_molecule.id in widget.molecules
        assert sample_molecule_methane.id in widget.molecules
        
        # 選択テスト
        widget.molecule_list.setCurrentRow(0)
        selected_id = widget.get_selected_molecule_id()
        selected_molecule = widget.get_selected_molecule()
        
        assert selected_id is not None
        assert selected_molecule is not None
        assert selected_molecule.name in [sample_molecule.name, sample_molecule_methane.name]
        
        # 削除テスト
        widget.remove_molecule(sample_molecule.id)
        assert widget.molecule_list.count() == 1
        assert sample_molecule.id not in widget.molecules
        
        widget.close()
    
    def test_simple_3d_viewer(self, qapp, sample_molecule):
        """Simple3DViewer のテスト"""
        from pyscf_front.gui.widgets.molecule_viewer import Simple3DViewer
        
        viewer = Simple3DViewer()
        
        # 初期状態の確認
        initial_text = viewer.viewer_area.toPlainText()
        assert "3D分子ビューアー" in initial_text
        assert "分子を読み込んでください" in initial_text
        
        # 分子を設定
        viewer.set_molecule(sample_molecule)
        
        # 表示内容の確認
        view_text = viewer.viewer_area.toPlainText()
        assert sample_molecule.name in view_text
        assert sample_molecule.get_molecular_formula() in view_text
        assert "重心:" in view_text
        assert "境界ボックス:" in view_text
        assert "原子一覧:" in view_text
        
        viewer.close()


class TestCalculationDialog:
    """CalculationDialog の GUI テスト"""
    
    def test_calculation_settings_dialog_initialization(self, qapp, sample_molecule, mock_pyscf):
        """CalculationSettingsDialog の初期化テスト"""
        with patch('pyscf_front.gui.dialogs.calculation_dialog.UnifiedCalculationEngine') as mock_engine_cls:
            mock_engine = Mock()
            mock_engine.get_available_methods.return_value = ['HF', 'B3LYP', 'PBE']
            mock_engine.get_available_basis_sets.return_value = ['STO-3G', '6-31G', '6-31G(d)']
            mock_engine.estimate_calculation_time.return_value = "約2分"
            
            from pyscf_front.gui.dialogs.calculation_dialog import CalculationSettingsDialog
            
            dialog = CalculationSettingsDialog(sample_molecule, mock_engine)
            
            # 基本的な属性の確認
            assert hasattr(dialog, 'method_combo')
            assert hasattr(dialog, 'basis_combo')
            assert hasattr(dialog, 'charge_spin')
            assert hasattr(dialog, 'mult_spin')
            assert hasattr(dialog, 'preview_text')
            
            # 初期値の確認
            assert dialog.charge_spin.value() == sample_molecule.charge
            assert dialog.mult_spin.value() == sample_molecule.multiplicity
            
            # コンボボックスの内容確認
            assert dialog.method_combo.count() > 0
            assert dialog.basis_combo.count() > 0
            
            dialog.close()
    
    def test_calculation_dialog_preview_update(self, qapp, qtbot, sample_molecule, mock_pyscf):
        """計算設定ダイアログのプレビュー更新テスト"""
        with patch('pyscf_front.gui.dialogs.calculation_dialog.UnifiedCalculationEngine') as mock_engine_cls:
            mock_engine = Mock()
            mock_engine.get_available_methods.return_value = ['HF', 'B3LYP']
            mock_engine.get_available_basis_sets.return_value = ['STO-3G', '6-31G']
            mock_engine.estimate_calculation_time.return_value = "約1分"
            
            from pyscf_front.gui.dialogs.calculation_dialog import CalculationSettingsDialog
            
            dialog = CalculationSettingsDialog(sample_molecule, mock_engine)
            qtbot.addWidget(dialog)
            
            # 設定を変更
            dialog.method_combo.setCurrentText('B3LYP')
            dialog.basis_combo.setCurrentText('6-31G')
            dialog.charge_spin.setValue(-1)
            dialog.mult_spin.setValue(2)
            
            # プレビューを更新
            dialog.update_preview()
            
            # プレビューテキストの確認
            preview_text = dialog.preview_text.toPlainText()
            assert "B3LYP" in preview_text
            assert "6-31G" in preview_text
            assert "電荷: -1" in preview_text
            assert "スピン多重度: 2" in preview_text
            assert sample_molecule.name in preview_text
            
            dialog.close()
    
    def test_calculation_start_confirmation(self, qapp, qtbot, sample_molecule, mock_pyscf):
        """計算開始確認ダイアログのテスト"""
        with patch('pyscf_front.gui.dialogs.calculation_dialog.CalculationEngine') as mock_engine_cls, \
             patch('pyscf_front.gui.dialogs.calculation_dialog.QMessageBox') as mock_msgbox:
            
            mock_engine = Mock()
            mock_engine.get_available_methods.return_value = ['HF']
            mock_engine.get_available_basis_sets.return_value = ['STO-3G']
            mock_engine.estimate_calculation_time.return_value = "約30秒"
            mock_engine.submit_calculation_with_persistence.return_value = "test-job-id"
            
            # メッセージボックスでYesを選択
            mock_msgbox.question.return_value = mock_msgbox.StandardButton.Yes
            mock_msgbox.information.return_value = None
            
            from pyscf_front.gui.dialogs.calculation_dialog import CalculationSettingsDialog
            
            dialog = CalculationSettingsDialog(sample_molecule, mock_engine)
            qtbot.addWidget(dialog)
            
            # start_calculationメソッドを直接呼び出し
            dialog.start_calculation()
            
            # 確認ダイアログが表示されたことを確認
            mock_msgbox.question.assert_called_once()
            
            # 計算が投入されたことを確認
            mock_engine.submit_calculation_with_persistence.assert_called_once()
            
            # 情報ダイアログが表示されたことを確認
            mock_msgbox.information.assert_called_once()
            
            dialog.close()


class TestGUIIntegration:
    """GUI 統合テスト"""
    
    def test_main_window_with_molecule_workflow(self, qapp, qtbot, sample_molecule, mock_pyscf):
        """メインウィンドウでの分子ワークフローテスト"""
        with patch('pyscf_front.gui.main_window.UnifiedCalculationEngine') as mock_engine_cls, \
             patch('pyscf_front.gui.main_window.MoleculeManager') as mock_mol_mgr_cls:
            
            mock_engine = Mock()
            mock_engine.resume_interrupted_calculations.return_value = []
            mock_engine.submit_calculation_with_persistence.return_value = "integration-test-job"
            mock_engine.list_jobs.return_value = []
            mock_engine.signals = Mock()
            mock_engine.signals.started = Mock()
            mock_engine.signals.progress = Mock()
            mock_engine.signals.completed = Mock()
            mock_engine.signals.failed = Mock()
            
            mock_engine_cls.return_value = mock_engine
            mock_mol_mgr_cls.return_value = Mock()
            
            from pyscf_front.gui.main_window import MainWindow
            
            window = MainWindow()
            qtbot.addWidget(window)
            
            # 1. 分子を設定
            window.molecule_widget.current_molecule = sample_molecule
            
            # 2. 分子選択イベントをシミュレート
            window.on_molecule_selected(sample_molecule)
            
            # 3. 計算設定が更新されたことを確認
            settings_text = window.calc_settings_text.toPlainText()
            assert sample_molecule.name in settings_text
            
            # 4. クイック計算実行
            qtbot.mouseClick(window.quick_calc_btn, Qt.MouseButton.LeftButton)
            
            # 5. 計算が投入されたことを確認
            mock_engine.submit_calculation_with_persistence.assert_called_once()
            
            # 6. ジョブリスト更新
            window.show_job_list()
            mock_engine.list_jobs.assert_called()
            
            window.close()
    
    def test_gui_error_handling(self, qapp, qtbot, mock_pyscf):
        """GUI エラーハンドリングテスト"""
        with patch('pyscf_front.gui.main_window.UnifiedCalculationEngine') as mock_engine_cls, \
             patch('pyscf_front.gui.main_window.MoleculeManager') as mock_mol_mgr_cls:
            
            mock_engine = Mock()
            mock_engine.resume_interrupted_calculations.return_value = []
            mock_engine.submit_calculation_with_persistence.side_effect = RuntimeError("計算エラー")
            mock_engine.signals = Mock()
            mock_engine.signals.started = Mock()
            mock_engine.signals.progress = Mock()
            mock_engine.signals.completed = Mock()
            mock_engine.signals.failed = Mock()
            
            mock_engine_cls.return_value = mock_engine
            mock_mol_mgr_cls.return_value = Mock()
            
            from pyscf_front.gui.main_window import MainWindow
            from pyscf_front.core.molecule import Molecule
            
            window = MainWindow()
            qtbot.addWidget(window)
            
            # 分子なしでクイック計算を試行
            window.molecule_widget.current_molecule = None
            qtbot.mouseClick(window.quick_calc_btn, Qt.MouseButton.LeftButton)
            
            # エラーメッセージがステータスバーに表示されることを確認
            status_message = window.status_bar.currentMessage()
            assert "分子を選択してください" in status_message
            
            # 分子を設定してエラーが発生する場合
            test_molecule = Molecule("Test")
            window.molecule_widget.current_molecule = test_molecule
            qtbot.mouseClick(window.quick_calc_btn, Qt.MouseButton.LeftButton)
            
            # エラーハンドリングが動作することを確認
            mock_engine.submit_calculation_with_persistence.assert_called()
            
            window.close()


if __name__ == "__main__":
    pytest.main([__file__, "-v", "--tb=short"])