"""
分子入力ダイアログ
"""
from typing import Optional
from PySide6.QtWidgets import (
    QDialog, QVBoxLayout, QHBoxLayout, QTabWidget, QWidget,
    QLabel, QLineEdit, QTextEdit, QPushButton, QSpinBox,
    QFormLayout, QMessageBox, QGroupBox,
    QFileDialog, QTableWidget, QTableWidgetItem, QHeaderView
)
from PySide6.QtCore import Signal, Slot
from PySide6.QtGui import QFont
from loguru import logger

from pyscf_front.core.molecule import Molecule, MoleculeBuilder, Atom


class MoleculeInputDialog(QDialog):
    """分子入力ダイアログ"""
    
    molecule_created = Signal(object)  # Molecule object
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("分子構造入力")
        self.setModal(True)
        self.setMinimumSize(600, 500)
        
        self.molecule = None
        self.setup_ui()
        self.connect_signals()
        
    def setup_ui(self):
        """UI設定"""
        layout = QVBoxLayout()
        
        # タブウィジェット
        self.tab_widget = QTabWidget()
        
        # SMILES入力タブ
        self.setup_smiles_tab()
        
        # XYZ入力タブ
        self.setup_xyz_tab()
        
        # 手動入力タブ
        self.setup_manual_tab()
        
        # プリセット分子タブ
        self.setup_preset_tab()
        
        layout.addWidget(self.tab_widget)
        
        # ボタン
        button_layout = QHBoxLayout()
        
        self.preview_btn = QPushButton("プレビュー")
        self.create_btn = QPushButton("作成")
        self.cancel_btn = QPushButton("キャンセル")
        
        self.create_btn.setDefault(True)
        
        button_layout.addStretch()
        button_layout.addWidget(self.preview_btn)
        button_layout.addWidget(self.create_btn)
        button_layout.addWidget(self.cancel_btn)
        
        layout.addLayout(button_layout)
        
        self.setLayout(layout)
    
    def setup_smiles_tab(self):
        """SMILESタブの設定"""
        widget = QWidget()
        layout = QVBoxLayout()
        
        # 説明
        desc_label = QLabel("SMILES記法を使用して分子構造を入力してください")
        desc_label.setWordWrap(True)
        layout.addWidget(desc_label)
        
        # 入力フォーム
        form_layout = QFormLayout()
        
        self.smiles_input = QLineEdit()
        self.smiles_input.setPlaceholderText("例: CCO (エタノール)")
        form_layout.addRow("SMILES:", self.smiles_input)
        
        self.smiles_name_input = QLineEdit()
        self.smiles_name_input.setPlaceholderText("分子名（オプション）")
        form_layout.addRow("分子名:", self.smiles_name_input)
        
        layout.addLayout(form_layout)
        
        # よく使用するSMILES例
        examples_group = QGroupBox("SMILES例")
        examples_layout = QVBoxLayout()
        
        examples = [
            ("水 (H2O)", "O"),
            ("メタン (CH4)", "C"),
            ("エタノール (C2H5OH)", "CCO"),
            ("ベンゼン (C6H6)", "c1ccccc1"),
            ("アセトン (C3H6O)", "CC(=O)C"),
            ("アンモニア (NH3)", "N")
        ]
        
        for name, smiles in examples:
            btn = QPushButton(f"{name}: {smiles}")
            btn.clicked.connect(lambda checked, s=smiles, n=name.split()[0]: self.set_smiles_example(s, n))
            examples_layout.addWidget(btn)
        
        examples_group.setLayout(examples_layout)
        layout.addWidget(examples_group)
        
        layout.addStretch()
        widget.setLayout(layout)
        self.tab_widget.addTab(widget, "SMILES")
    
    def setup_xyz_tab(self):
        """XYZタブの設定"""
        widget = QWidget()
        layout = QVBoxLayout()
        
        # 説明
        desc_label = QLabel("XYZ形式で分子構造を入力してください")
        layout.addWidget(desc_label)
        
        # 入力フォーム
        form_layout = QFormLayout()
        
        self.xyz_name_input = QLineEdit()
        self.xyz_name_input.setPlaceholderText("分子名")
        form_layout.addRow("分子名:", self.xyz_name_input)
        
        layout.addLayout(form_layout)
        
        # XYZ入力エリア
        self.xyz_text = QTextEdit()
        self.xyz_text.setPlaceholderText("""3
Water molecule
O    0.000000    0.000000    0.000000
H    0.757000    0.586000    0.000000
H   -0.757000    0.586000    0.000000""")
        self.xyz_text.setFont(QFont("Courier", 10))
        layout.addWidget(self.xyz_text)
        
        # ファイル読み込みボタン
        file_layout = QHBoxLayout()
        self.load_xyz_btn = QPushButton("XYZファイル読み込み")
        file_layout.addWidget(self.load_xyz_btn)
        file_layout.addStretch()
        layout.addLayout(file_layout)
        
        widget.setLayout(layout)
        self.tab_widget.addTab(widget, "XYZ")
    
    def setup_manual_tab(self):
        """手動入力タブの設定"""
        widget = QWidget()
        layout = QVBoxLayout()
        
        # 基本情報
        info_group = QGroupBox("分子情報")
        info_layout = QFormLayout()
        
        self.manual_name_input = QLineEdit()
        self.manual_name_input.setPlaceholderText("分子名")
        info_layout.addRow("分子名:", self.manual_name_input)
        
        self.charge_spin = QSpinBox()
        self.charge_spin.setRange(-10, 10)
        self.charge_spin.setValue(0)
        info_layout.addRow("電荷:", self.charge_spin)
        
        self.multiplicity_spin = QSpinBox()
        self.multiplicity_spin.setRange(1, 10)
        self.multiplicity_spin.setValue(1)
        info_layout.addRow("スピン多重度:", self.multiplicity_spin)
        
        info_group.setLayout(info_layout)
        layout.addWidget(info_group)
        
        # 原子入力テーブル
        atoms_group = QGroupBox("原子")
        atoms_layout = QVBoxLayout()
        
        self.atoms_table = QTableWidget()
        self.atoms_table.setColumnCount(4)
        self.atoms_table.setHorizontalHeaderLabels(["元素", "X座標", "Y座標", "Z座標"])
        self.atoms_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        
        atoms_layout.addWidget(self.atoms_table)
        
        # 原子操作ボタン
        atoms_btn_layout = QHBoxLayout()
        self.add_atom_btn = QPushButton("原子追加")
        self.remove_atom_btn = QPushButton("原子削除")
        
        atoms_btn_layout.addWidget(self.add_atom_btn)
        atoms_btn_layout.addWidget(self.remove_atom_btn)
        atoms_btn_layout.addStretch()
        
        atoms_layout.addLayout(atoms_btn_layout)
        atoms_group.setLayout(atoms_layout)
        layout.addWidget(atoms_group)
        
        widget.setLayout(layout)
        self.tab_widget.addTab(widget, "手動入力")
    
    def setup_preset_tab(self):
        """プリセット分子タブの設定"""
        widget = QWidget()
        layout = QVBoxLayout()
        
        desc_label = QLabel("よく使用される分子構造から選択してください")
        layout.addWidget(desc_label)
        
        # プリセット分子リスト
        presets = [
            ("水 (H2O)", "water"),
            ("メタン (CH4)", "methane"),
            ("ベンゼン (C6H6)", "benzene"),
            ("アンモニア (NH3)", "ammonia"),
            ("二酸化炭素 (CO2)", "co2"),
            ("エタン (C2H6)", "ethane")
        ]
        
        for name, preset_id in presets:
            btn = QPushButton(name)
            btn.clicked.connect(lambda checked, pid=preset_id: self.select_preset(pid))
            layout.addWidget(btn)
        
        layout.addStretch()
        widget.setLayout(layout)
        self.tab_widget.addTab(widget, "プリセット")
    
    def connect_signals(self):
        """シグナル接続"""
        self.preview_btn.clicked.connect(self.preview_molecule)
        self.create_btn.clicked.connect(self.create_molecule)
        self.cancel_btn.clicked.connect(self.reject)
        
        self.load_xyz_btn.clicked.connect(self.load_xyz_file)
        self.add_atom_btn.clicked.connect(self.add_atom_row)
        self.remove_atom_btn.clicked.connect(self.remove_atom_row)
    
    @Slot(str, str)
    def set_smiles_example(self, smiles: str, name: str):
        """SMILES例を設定"""
        self.smiles_input.setText(smiles)
        self.smiles_name_input.setText(name)
    
    @Slot(str)
    def select_preset(self, preset_id: str):
        """プリセット分子を選択"""
        try:
            if preset_id == "water":
                self.molecule = MoleculeBuilder.create_water()
            elif preset_id == "methane":
                self.molecule = MoleculeBuilder.create_methane()
            elif preset_id == "benzene":
                self.molecule = MoleculeBuilder.create_benzene()
            else:
                # その他のプリセットはSMILESで作成
                smiles_map = {
                    "ammonia": "N",
                    "co2": "O=C=O",
                    "ethane": "CC"
                }
                if preset_id in smiles_map:
                    self.molecule = MoleculeBuilder.from_smiles(smiles_map[preset_id])
            
            if self.molecule:
                QMessageBox.information(self, "成功", f"プリセット分子 '{self.molecule.name}' を選択しました")
                logger.info(f"Selected preset molecule: {self.molecule.name}")
            else:
                QMessageBox.warning(self, "エラー", "プリセット分子の作成に失敗しました")
                
        except Exception as e:
            QMessageBox.critical(self, "エラー", f"プリセット分子の作成中にエラーが発生しました: {e}")
            logger.error(f"Error creating preset molecule {preset_id}: {e}")
    
    @Slot()
    def load_xyz_file(self):
        """XYZファイルを読み込み"""
        file_path, _ = QFileDialog.getOpenFileName(
            self, "XYZファイルを開く", "", "XYZ Files (*.xyz);;All Files (*)"
        )
        
        if file_path:
            try:
                with open(file_path, 'r', encoding='utf-8') as f:
                    content = f.read()
                self.xyz_text.setPlainText(content)
                
                # ファイル名から分子名を推測
                import os
                filename = os.path.splitext(os.path.basename(file_path))[0]
                self.xyz_name_input.setText(filename)
                
                logger.info(f"Loaded XYZ file: {file_path}")
                
            except Exception as e:
                QMessageBox.critical(self, "エラー", f"ファイルの読み込みに失敗しました: {e}")
                logger.error(f"Error loading XYZ file {file_path}: {e}")
    
    @Slot()
    def add_atom_row(self):
        """原子行を追加"""
        row = self.atoms_table.rowCount()
        self.atoms_table.insertRow(row)
        
        # デフォルト値を設定
        self.atoms_table.setItem(row, 0, QTableWidgetItem("C"))
        self.atoms_table.setItem(row, 1, QTableWidgetItem("0.0"))
        self.atoms_table.setItem(row, 2, QTableWidgetItem("0.0"))
        self.atoms_table.setItem(row, 3, QTableWidgetItem("0.0"))
    
    @Slot()
    def remove_atom_row(self):
        """原子行を削除"""
        current_row = self.atoms_table.currentRow()
        if current_row >= 0:
            self.atoms_table.removeRow(current_row)
    
    @Slot()
    def preview_molecule(self):
        """分子をプレビュー"""
        molecule = self._create_molecule_from_current_tab()
        if molecule:
            info = f"""
分子名: {molecule.name}
分子式: {molecule.get_molecular_formula()}
原子数: {len(molecule.atoms)}
電荷: {molecule.charge}
スピン多重度: {molecule.multiplicity}

原子座標:
{molecule.to_xyz()}
            """
            QMessageBox.information(self, "分子プレビュー", info)
        else:
            QMessageBox.warning(self, "エラー", "分子の作成に失敗しました")
    
    @Slot()
    def create_molecule(self):
        """分子を作成"""
        molecule = self._create_molecule_from_current_tab()
        if molecule:
            self.molecule = molecule
            self.molecule_created.emit(self.molecule)
            self.accept()
        else:
            QMessageBox.warning(self, "エラー", "分子の作成に失敗しました")
    
    def _create_molecule_from_current_tab(self) -> Optional[Molecule]:
        """現在のタブから分子を作成"""
        current_tab = self.tab_widget.currentIndex()
        
        try:
            if current_tab == 0:  # SMILES
                return self._create_from_smiles()
            elif current_tab == 1:  # XYZ
                return self._create_from_xyz()
            elif current_tab == 2:  # Manual
                return self._create_from_manual()
            elif current_tab == 3:  # Preset
                return self.molecule
            
        except Exception as e:
            logger.error(f"Error creating molecule from tab {current_tab}: {e}")
            QMessageBox.critical(self, "エラー", f"分子作成中にエラーが発生しました: {e}")
        
        return None
    
    def _create_from_smiles(self) -> Optional[Molecule]:
        """SMILES入力から分子を作成"""
        smiles = self.smiles_input.text().strip()
        name = self.smiles_name_input.text().strip()
        
        if not smiles:
            QMessageBox.warning(self, "入力エラー", "SMILES文字列を入力してください")
            return None
        
        molecule = MoleculeBuilder.from_smiles(smiles, name)
        if not molecule:
            QMessageBox.warning(self, "エラー", "無効なSMILES文字列です")
        
        return molecule
    
    def _create_from_xyz(self) -> Optional[Molecule]:
        """XYZ入力から分子を作成"""
        xyz_text = self.xyz_text.toPlainText().strip()
        name = self.xyz_name_input.text().strip()
        
        if not xyz_text:
            QMessageBox.warning(self, "入力エラー", "XYZ座標を入力してください")
            return None
        
        molecule = Molecule(name or "XYZ_Molecule")
        molecule.from_xyz_string(xyz_text)
        return molecule
    
    def _create_from_manual(self) -> Optional[Molecule]:
        """手動入力から分子を作成"""
        name = self.manual_name_input.text().strip() or "Manual_Molecule"
        charge = self.charge_spin.value()
        multiplicity = self.multiplicity_spin.value()
        
        if self.atoms_table.rowCount() == 0:
            QMessageBox.warning(self, "入力エラー", "少なくとも1つの原子を追加してください")
            return None
        
        molecule = Molecule(name, charge, multiplicity)
        
        for row in range(self.atoms_table.rowCount()):
            try:
                symbol_item = self.atoms_table.item(row, 0)
                x_item = self.atoms_table.item(row, 1)
                y_item = self.atoms_table.item(row, 2)
                z_item = self.atoms_table.item(row, 3)
                
                if not all([symbol_item, x_item, y_item, z_item]):
                    QMessageBox.warning(self, "入力エラー", f"行 {row + 1} のデータが不完全です")
                    return None
                
                # Type assertion after null check
                assert symbol_item is not None
                assert x_item is not None
                assert y_item is not None
                assert z_item is not None
                
                symbol = symbol_item.text().strip()
                x = float(x_item.text())
                y = float(y_item.text())
                z = float(z_item.text())
                
                molecule.add_atom(Atom(symbol, x, y, z))
                
            except (ValueError, AttributeError) as e:
                QMessageBox.warning(self, "入力エラー", f"行 {row + 1} のデータが無効です: {e}")
                return None
        
        return molecule