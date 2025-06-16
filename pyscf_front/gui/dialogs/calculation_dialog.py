"""
計算設定ダイアログ
"""
from PySide6.QtWidgets import (
    QDialog, QVBoxLayout, QHBoxLayout, QFormLayout, QGroupBox,
    QLabel, QComboBox, QSpinBox, QTextEdit, QPushButton,
    QMessageBox, QProgressBar, QCheckBox
)
from PySide6.QtCore import Signal, Slot
from PySide6.QtGui import QFont
from loguru import logger

from pyscf_front.core.molecule import Molecule
from pyscf_front.core.calculation_engine_unified import UnifiedCalculationEngine


class CalculationSettingsDialog(QDialog):
    """計算設定ダイアログ"""
    
    calculation_requested = Signal(str, str, str, int, int, dict)  # molecule_id, method, basis, charge, mult, params
    
    def __init__(self, molecule: Molecule, calculation_engine: "UnifiedCalculationEngine", parent=None):  # type: ignore
        super().__init__(parent)
        self.molecule = molecule
        self.calculation_engine = calculation_engine
        
        self.setWindowTitle(f"計算設定 - {molecule.name}")
        self.setModal(True)
        self.setMinimumSize(500, 600)
        
        self.setup_ui()
        self.connect_signals()
        self.load_defaults()
        
    def setup_ui(self):
        """UI設定"""
        layout = QVBoxLayout()
        
        # 分子情報
        mol_group = QGroupBox("分子情報")
        mol_layout = QFormLayout()
        
        mol_layout.addRow("分子名:", QLabel(self.molecule.name))
        mol_layout.addRow("分子式:", QLabel(self.molecule.get_molecular_formula()))
        mol_layout.addRow("原子数:", QLabel(str(len(self.molecule.atoms))))
        
        mol_group.setLayout(mol_layout)
        layout.addWidget(mol_group)
        
        # 計算設定
        calc_group = QGroupBox("計算設定")
        calc_layout = QFormLayout()
        
        # 計算手法
        self.method_combo = QComboBox()
        self.method_combo.addItems(self.calculation_engine.get_available_methods())
        calc_layout.addRow("計算手法:", self.method_combo)
        
        # 基底関数
        self.basis_combo = QComboBox()
        self.basis_combo.addItems(self.calculation_engine.get_available_basis_sets())
        calc_layout.addRow("基底関数:", self.basis_combo)
        
        # 電荷
        self.charge_spin = QSpinBox()
        self.charge_spin.setRange(-10, 10)
        self.charge_spin.setValue(self.molecule.charge)
        calc_layout.addRow("電荷:", self.charge_spin)
        
        # スピン多重度
        self.mult_spin = QSpinBox()
        self.mult_spin.setRange(1, 10)
        self.mult_spin.setValue(self.molecule.multiplicity)
        calc_layout.addRow("スピン多重度:", self.mult_spin)
        
        calc_group.setLayout(calc_layout)
        layout.addWidget(calc_group)
        
        # 詳細設定
        advanced_group = QGroupBox("詳細設定")
        advanced_layout = QVBoxLayout()
        
        # 収束設定
        convergence_layout = QFormLayout()
        
        self.max_cycles_spin = QSpinBox()
        self.max_cycles_spin.setRange(10, 1000)
        self.max_cycles_spin.setValue(50)
        convergence_layout.addRow("最大サイクル数:", self.max_cycles_spin)
        
        self.conv_threshold_combo = QComboBox()
        self.conv_threshold_combo.addItems(["1e-6", "1e-8", "1e-10", "1e-12"])
        self.conv_threshold_combo.setCurrentText("1e-8")
        convergence_layout.addRow("収束判定値:", self.conv_threshold_combo)
        
        advanced_layout.addLayout(convergence_layout)
        
        # その他のオプション
        options_layout = QVBoxLayout()
        
        self.verbose_check = QCheckBox("詳細出力")
        self.save_orbitals_check = QCheckBox("軌道を保存")
        
        options_layout.addWidget(self.verbose_check)
        options_layout.addWidget(self.save_orbitals_check)
        
        advanced_layout.addLayout(options_layout)
        advanced_group.setLayout(advanced_layout)
        layout.addWidget(advanced_group)
        
        # 推定時間表示
        self.time_label = QLabel("推定計算時間: -")
        self.time_label.setStyleSheet("font-weight: bold; color: #0d7377;")
        layout.addWidget(self.time_label)
        
        # プレビューエリア
        preview_group = QGroupBox("設定プレビュー")
        preview_layout = QVBoxLayout()
        
        self.preview_text = QTextEdit()
        self.preview_text.setReadOnly(True)
        self.preview_text.setFont(QFont("Courier", 9))
        self.preview_text.setMaximumHeight(150)
        
        preview_layout.addWidget(self.preview_text)
        preview_group.setLayout(preview_layout)
        layout.addWidget(preview_group)
        
        # ボタン
        button_layout = QHBoxLayout()
        
        self.estimate_btn = QPushButton("時間推定")
        self.start_btn = QPushButton("計算開始")
        self.cancel_btn = QPushButton("キャンセル")
        
        self.start_btn.setDefault(True)
        self.start_btn.setStyleSheet("background-color: #0d7377; color: white; font-weight: bold;")
        
        button_layout.addWidget(self.estimate_btn)
        button_layout.addStretch()
        button_layout.addWidget(self.start_btn)
        button_layout.addWidget(self.cancel_btn)
        
        layout.addLayout(button_layout)
        
        self.setLayout(layout)
    
    def connect_signals(self):
        """シグナル接続"""
        self.method_combo.currentTextChanged.connect(self.update_preview)
        self.basis_combo.currentTextChanged.connect(self.update_preview)
        self.charge_spin.valueChanged.connect(self.update_preview)
        self.mult_spin.valueChanged.connect(self.update_preview)
        
        self.estimate_btn.clicked.connect(self.estimate_time)
        self.start_btn.clicked.connect(self.start_calculation)
        self.cancel_btn.clicked.connect(self.reject)
    
    def load_defaults(self):
        """デフォルト値を読み込み"""
        # 分子サイズに応じた推奨設定
        atom_count = len(self.molecule.atoms)
        
        if atom_count <= 10:
            # 小さい分子: 高精度設定
            self.method_combo.setCurrentText("B3LYP")
            self.basis_combo.setCurrentText("6-31G(d,p)")
        elif atom_count <= 30:
            # 中程度の分子: バランス設定
            self.method_combo.setCurrentText("B3LYP")
            self.basis_combo.setCurrentText("6-31G(d)")
        else:
            # 大きい分子: 効率重視設定
            self.method_combo.setCurrentText("PBE")
            self.basis_combo.setCurrentText("6-31G")
        
        self.update_preview()
        self.estimate_time()
    
    @Slot()
    def update_preview(self):
        """プレビューを更新"""
        method = self.method_combo.currentText()
        basis = self.basis_combo.currentText()
        charge = self.charge_spin.value()
        mult = self.mult_spin.value()
        
        preview_text = f"""計算設定プレビュー:

分子: {self.molecule.name}
分子式: {self.molecule.get_molecular_formula()}
原子数: {len(self.molecule.atoms)}

計算設定:
  手法: {method}
  基底関数: {basis}
  電荷: {charge}
  スピン多重度: {mult}
  
詳細設定:
  最大サイクル数: {self.max_cycles_spin.value()}
  収束判定値: {self.conv_threshold_combo.currentText()}
  詳細出力: {'はい' if self.verbose_check.isChecked() else 'いいえ'}
  軌道保存: {'はい' if self.save_orbitals_check.isChecked() else 'いいえ'}

座標情報:
{self.molecule.to_xyz()}
"""
        
        self.preview_text.setPlainText(preview_text)
    
    @Slot()
    def estimate_time(self):
        """計算時間を推定"""
        method = self.method_combo.currentText()
        basis = self.basis_combo.currentText()
        
        estimated_time = self.calculation_engine.estimate_calculation_time(
            self.molecule, method, basis
        )
        
        self.time_label.setText(f"推定計算時間: {estimated_time}")
        logger.info(f"Estimated calculation time: {estimated_time}")
    
    @Slot()
    def start_calculation(self):
        """計算を開始"""
        method = self.method_combo.currentText()
        basis = self.basis_combo.currentText()
        charge = self.charge_spin.value()
        mult = self.mult_spin.value()
        
        # パラメータ辞書作成
        parameters = {
            'max_cycles': self.max_cycles_spin.value(),
            'conv_threshold': float(self.conv_threshold_combo.currentText()),
            'verbose': self.verbose_check.isChecked(),
            'save_orbitals': self.save_orbitals_check.isChecked()
        }
        
        # 確認ダイアログ
        reply = QMessageBox.question(
            self,
            "計算開始確認",
            f"""以下の設定で計算を開始しますか？

分子: {self.molecule.name}
手法: {method}
基底関数: {basis}
推定時間: {self.time_label.text().split(': ')[1]}

計算は別スレッドで実行されるため、アプリケーションの操作は継続できます。""",
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
            QMessageBox.StandardButton.Yes
        )
        
        if reply == QMessageBox.StandardButton.Yes:
            try:
                # 計算エンジンに投入（データベース永続化付き）
                if hasattr(self.calculation_engine, 'submit_calculation_with_persistence'):
                    # 統合版の場合
                    job_id = self.calculation_engine.submit_calculation_with_persistence(
                        self.molecule,
                        method=method,
                        basis_set=basis,
                        charge=charge,
                        multiplicity=mult,
                        parameters=parameters,
                        instance_name=f"Calc_{self.molecule.name}_{method}_{basis}"
                    )
                else:
                    # 従来版の場合
                    job_id = self.calculation_engine.submit_calculation(
                        self.molecule,
                        method=method,
                        basis_set=basis,
                        charge=charge,
                        multiplicity=mult,
                        parameters=parameters
                    )
                
                QMessageBox.information(
                    self,
                    "計算開始",
                    f"計算を開始しました。\nジョブID: {job_id}\n\n進捗はメインウィンドウで確認できます。"
                )
                
                self.accept()
                
            except Exception as e:
                QMessageBox.critical(
                    self,
                    "エラー",
                    f"計算の開始に失敗しました:\n{str(e)}"
                )
                logger.error(f"Failed to start calculation: {e}")


class CalculationProgressDialog(QDialog):
    """計算進捗ダイアログ"""
    
    def __init__(self, job_id: str, molecule_name: str, parent=None):
        super().__init__(parent)
        self.job_id = job_id
        
        self.setWindowTitle(f"計算進捗 - {molecule_name}")
        self.setModal(False)
        self.setMinimumSize(400, 200)
        
        self.setup_ui()
    
    def setup_ui(self):
        """UI設定"""
        layout = QVBoxLayout()
        
        # 情報表示
        self.info_label = QLabel(f"ジョブID: {self.job_id}")
        layout.addWidget(self.info_label)
        
        # 進捗バー
        self.progress_bar = QProgressBar()
        self.progress_bar.setRange(0, 100)
        self.progress_bar.setValue(0)
        layout.addWidget(self.progress_bar)
        
        # メッセージ
        self.message_label = QLabel("計算を準備中...")
        layout.addWidget(self.message_label)
        
        # ボタン
        button_layout = QHBoxLayout()
        
        self.cancel_btn = QPushButton("キャンセル")
        self.close_btn = QPushButton("閉じる")
        
        button_layout.addStretch()
        button_layout.addWidget(self.cancel_btn)
        button_layout.addWidget(self.close_btn)
        
        layout.addLayout(button_layout)
        
        # 接続
        self.cancel_btn.clicked.connect(self.cancel_calculation)
        self.close_btn.clicked.connect(self.close)
        
        self.setLayout(layout)
    
    @Slot(float, str)
    def update_progress(self, progress: float, message: str):
        """進捗を更新"""
        self.progress_bar.setValue(int(progress * 100))
        self.message_label.setText(message)
    
    @Slot()
    def calculation_completed(self):
        """計算完了時の処理"""
        self.progress_bar.setValue(100)
        self.message_label.setText("計算完了!")
        self.cancel_btn.setEnabled(False)
    
    @Slot()
    def calculation_failed(self, error_msg: str):
        """計算失敗時の処理"""
        self.message_label.setText(f"計算失敗: {error_msg}")
        self.cancel_btn.setEnabled(False)
    
    @Slot()
    def cancel_calculation(self):
        """計算をキャンセル"""
        # 実装は親ウィンドウで処理
        self.reject()