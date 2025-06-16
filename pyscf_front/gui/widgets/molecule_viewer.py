"""
分子表示ウィジェット
"""
from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton,
    QTextEdit, QGroupBox, QListWidget, QListWidgetItem
)
from PySide6.QtCore import Qt, Signal, Slot
from PySide6.QtGui import QFont
from loguru import logger
from typing import Optional, Dict

from pyscf_front.core.molecule import Molecule


class MoleculeInfoWidget(QWidget):
    """分子情報表示ウィジェット"""
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.molecule = None
        self.setup_ui()
    
    def setup_ui(self):
        """UI設定"""
        layout = QVBoxLayout()
        
        # 基本情報グループ
        info_group = QGroupBox("分子情報")
        info_layout = QVBoxLayout()
        
        self.name_label = QLabel("名前: -")
        self.formula_label = QLabel("分子式: -")
        self.atoms_label = QLabel("原子数: -")
        self.charge_label = QLabel("電荷: -")
        self.multiplicity_label = QLabel("スピン多重度: -")
        
        info_layout.addWidget(self.name_label)
        info_layout.addWidget(self.formula_label)
        info_layout.addWidget(self.atoms_label)
        info_layout.addWidget(self.charge_label)
        info_layout.addWidget(self.multiplicity_label)
        
        info_group.setLayout(info_layout)
        layout.addWidget(info_group)
        
        # 座標表示グループ
        coords_group = QGroupBox("原子座標")
        coords_layout = QVBoxLayout()
        
        self.coords_text = QTextEdit()
        self.coords_text.setReadOnly(True)
        self.coords_text.setFont(QFont("Courier", 9))
        self.coords_text.setMaximumHeight(200)
        
        coords_layout.addWidget(self.coords_text)
        coords_group.setLayout(coords_layout)
        layout.addWidget(coords_group)
        
        layout.addStretch()
        self.setLayout(layout)
    
    def set_molecule(self, molecule: Molecule):
        """分子を設定"""
        self.molecule = molecule
        self.update_display()
    
    def update_display(self):
        """表示を更新"""
        if self.molecule:
            self.name_label.setText(f"名前: {self.molecule.name}")
            self.formula_label.setText(f"分子式: {self.molecule.get_molecular_formula()}")
            self.atoms_label.setText(f"原子数: {len(self.molecule.atoms)}")
            self.charge_label.setText(f"電荷: {self.molecule.charge}")
            self.multiplicity_label.setText(f"スピン多重度: {self.molecule.multiplicity}")
            
            # 座標表示
            coords_text = "原子    X座標       Y座標       Z座標\n"
            coords_text += "-" * 40 + "\n"
            for atom in self.molecule.atoms:
                coords_text += f"{atom.symbol:2s} {atom.x:10.6f} {atom.y:10.6f} {atom.z:10.6f}\n"
            
            self.coords_text.setPlainText(coords_text)
        else:
            self.clear_display()
    
    def clear_display(self):
        """表示をクリア"""
        self.name_label.setText("名前: -")
        self.formula_label.setText("分子式: -")
        self.atoms_label.setText("原子数: -")
        self.charge_label.setText("電荷: -")
        self.multiplicity_label.setText("スピン多重度: -")
        self.coords_text.clear()


class Simple3DViewer(QWidget):
    """簡易3D分子ビューアー（テキストベース）"""
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.molecule = None
        self.setup_ui()
    
    def setup_ui(self):
        """UI設定"""
        layout = QVBoxLayout()
        
        # ヘッダー
        header_layout = QHBoxLayout()
        title_label = QLabel("3D分子構造")
        title_label.setStyleSheet("font-weight: bold; font-size: 14px;")
        
        self.view_controls = QHBoxLayout()
        self.center_btn = QPushButton("中央表示")
        self.reset_btn = QPushButton("リセット")
        
        self.view_controls.addWidget(self.center_btn)
        self.view_controls.addWidget(self.reset_btn)
        self.view_controls.addStretch()
        
        header_layout.addWidget(title_label)
        header_layout.addStretch()
        header_layout.addLayout(self.view_controls)
        
        layout.addLayout(header_layout)
        
        # 3D表示エリア（プレースホルダー）
        self.viewer_area = QTextEdit()
        self.viewer_area.setReadOnly(True)
        self.viewer_area.setFont(QFont("Courier", 10))
        self.viewer_area.setPlainText("3D分子ビューアー\n\n分子を読み込んでください...")
        
        layout.addWidget(self.viewer_area)
        
        # 接続
        self.center_btn.clicked.connect(self.center_view)
        self.reset_btn.clicked.connect(self.reset_view)
        
        self.setLayout(layout)
    
    def set_molecule(self, molecule: Optional[Molecule]):
        """分子を設定"""
        self.molecule = molecule
        self.update_view()
    
    def update_view(self):
        """ビューを更新"""
        if not self.molecule:
            self.viewer_area.setPlainText("3D分子ビューアー\n\n分子を読み込んでください...")
            return
        
        # 簡易的なテキスト表示
        view_text = f"分子: {self.molecule.name}\n"
        view_text += f"分子式: {self.molecule.get_molecular_formula()}\n\n"
        
        # 重心計算
        center = self.molecule.calculate_center_of_mass()
        view_text += f"重心: ({center[0]:.3f}, {center[1]:.3f}, {center[2]:.3f})\n\n"
        
        # 境界ボックス
        min_coords, max_coords = self.molecule.get_bounds()
        view_text += f"境界ボックス:\n"
        view_text += f"  最小: ({min_coords[0]:.3f}, {min_coords[1]:.3f}, {min_coords[2]:.3f})\n"
        view_text += f"  最大: ({max_coords[0]:.3f}, {max_coords[1]:.3f}, {max_coords[2]:.3f})\n\n"
        
        # 原子リスト
        view_text += "原子一覧:\n"
        view_text += "-" * 50 + "\n"
        for i, atom in enumerate(self.molecule.atoms):
            distance_from_center = ((atom.x - center[0])**2 + 
                                  (atom.y - center[1])**2 + 
                                  (atom.z - center[2])**2)**0.5
            view_text += f"{i+1:3d}. {atom.symbol:2s} "
            view_text += f"({atom.x:7.3f}, {atom.y:7.3f}, {atom.z:7.3f}) "
            view_text += f"距離: {distance_from_center:.3f}\n"
        
        # 結合情報
        if self.molecule.bonds:
            view_text += "\n結合情報:\n"
            view_text += "-" * 30 + "\n"
            for i, (atom1, atom2, order) in enumerate(self.molecule.bonds):
                atom1_symbol = self.molecule.atoms[atom1].symbol
                atom2_symbol = self.molecule.atoms[atom2].symbol
                view_text += f"{i+1:3d}. {atom1_symbol}{atom1+1}-{atom2_symbol}{atom2+1} (結合次数: {order})\n"
        
        self.viewer_area.setPlainText(view_text)
    
    @Slot()
    def center_view(self):
        """ビューを中央に"""
        logger.info("Center view requested")
        if self.molecule:
            self.update_view()
    
    @Slot()
    def reset_view(self):
        """ビューをリセット"""
        logger.info("Reset view requested")
        if self.molecule:
            self.update_view()


class MoleculeListWidget(QWidget):
    """分子リスト表示ウィジェット"""
    
    molecule_selected = Signal(str)  # molecule_id
    molecule_double_clicked = Signal(str)  # molecule_id
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.molecules: Dict[str, Molecule] = {}
        self.setup_ui()
        self.connect_signals()
    
    def setup_ui(self):
        """UI設定"""
        layout = QVBoxLayout()
        
        # ヘッダー
        header_layout = QHBoxLayout()
        title_label = QLabel("分子リスト")
        title_label.setStyleSheet("font-weight: bold; font-size: 14px;")
        
        self.add_btn = QPushButton("追加")
        self.remove_btn = QPushButton("削除")
        
        header_layout.addWidget(title_label)
        header_layout.addStretch()
        header_layout.addWidget(self.add_btn)
        header_layout.addWidget(self.remove_btn)
        
        layout.addLayout(header_layout)
        
        # 分子リスト
        self.molecule_list = QListWidget()
        layout.addWidget(self.molecule_list)
        
        self.setLayout(layout)
    
    def connect_signals(self):
        """シグナル接続"""
        self.molecule_list.itemSelectionChanged.connect(self.on_selection_changed)
        self.molecule_list.itemDoubleClicked.connect(self.on_item_double_clicked)
    
    def add_molecule(self, molecule: Molecule):
        """分子を追加"""
        self.molecules[molecule.id] = molecule
        
        # リストアイテム作成
        item_text = f"{molecule.name} ({molecule.get_molecular_formula()})"
        item = QListWidgetItem(item_text)
        item.setData(Qt.ItemDataRole.UserRole, molecule.id)
        
        self.molecule_list.addItem(item)
        logger.info(f"Added molecule to list: {molecule.name}")
    
    def remove_molecule(self, molecule_id: str):
        """分子を削除"""
        if molecule_id in self.molecules:
            del self.molecules[molecule_id]
            
            # リストアイテムを削除
            for i in range(self.molecule_list.count()):
                item = self.molecule_list.item(i)
                if item.data(Qt.ItemDataRole.UserRole) == molecule_id:
                    self.molecule_list.takeItem(i)
                    break
            
            logger.info(f"Removed molecule from list: {molecule_id}")
    
    def get_selected_molecule_id(self) -> Optional[str]:
        """選択された分子IDを取得"""
        current_item = self.molecule_list.currentItem()
        if current_item:
            return current_item.data(Qt.ItemDataRole.UserRole)
        return None
    
    def get_selected_molecule(self) -> Optional[Molecule]:
        """選択された分子を取得"""
        molecule_id = self.get_selected_molecule_id()
        if molecule_id:
            return self.molecules.get(molecule_id)
        return None
    
    @Slot()
    def on_selection_changed(self):
        """選択変更時の処理"""
        molecule_id = self.get_selected_molecule_id()
        if molecule_id:
            self.molecule_selected.emit(molecule_id)
    
    @Slot(QListWidgetItem)
    def on_item_double_clicked(self, item: QListWidgetItem):
        """アイテムダブルクリック時の処理"""
        molecule_id = item.data(Qt.ItemDataRole.UserRole)
        if molecule_id:
            self.molecule_double_clicked.emit(molecule_id)


class MoleculeManagerWidget(QWidget):
    """分子管理統合ウィジェット"""
    
    molecule_selected = Signal(object)  # Molecule
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.current_molecule = None
        self.setup_ui()
        self.connect_signals()
    
    def setup_ui(self):
        """UI設定"""
        layout = QHBoxLayout()
        
        # 左側: 分子リスト
        left_widget = QWidget()
        left_layout = QVBoxLayout()
        
        self.molecule_list = MoleculeListWidget()
        left_layout.addWidget(self.molecule_list)
        
        left_widget.setLayout(left_layout)
        left_widget.setMaximumWidth(300)
        
        # 右側: 分子情報と3Dビューアー
        right_widget = QWidget()
        right_layout = QVBoxLayout()
        
        # 分子情報
        self.molecule_info = MoleculeInfoWidget()
        right_layout.addWidget(self.molecule_info)
        
        # 3Dビューアー
        self.viewer_3d = Simple3DViewer()
        right_layout.addWidget(self.viewer_3d)
        
        right_widget.setLayout(right_layout)
        
        # レイアウトに追加
        layout.addWidget(left_widget)
        layout.addWidget(right_widget)
        
        self.setLayout(layout)
    
    def connect_signals(self):
        """シグナル接続"""
        self.molecule_list.molecule_selected.connect(self.on_molecule_selected)
        self.molecule_list.add_btn.clicked.connect(self.add_molecule_dialog)
        self.molecule_list.remove_btn.clicked.connect(self.remove_selected_molecule)
    
    def add_molecule(self, molecule: Molecule):
        """分子を追加"""
        self.molecule_list.add_molecule(molecule)
    
    @Slot(str)
    def on_molecule_selected(self, molecule_id: str):
        """分子選択時の処理"""
        molecule = self.molecule_list.molecules.get(molecule_id)
        if molecule:
            self.current_molecule = molecule
            self.molecule_info.set_molecule(molecule)
            self.viewer_3d.set_molecule(molecule)
            self.molecule_selected.emit(molecule)
            logger.info(f"Selected molecule: {molecule.name}")
    
    @Slot()
    def add_molecule_dialog(self):
        """分子追加ダイアログを表示"""
        from pyscf_front.gui.dialogs.molecule_input_dialog import MoleculeInputDialog
        
        dialog = MoleculeInputDialog(self)
        dialog.molecule_created.connect(self.add_molecule)
        dialog.exec()
    
    @Slot()
    def remove_selected_molecule(self):
        """選択された分子を削除"""
        molecule_id = self.molecule_list.get_selected_molecule_id()
        if molecule_id:
            self.molecule_list.remove_molecule(molecule_id)
            self.molecule_info.clear_display()
            self.viewer_3d.set_molecule(None)
            self.current_molecule = None