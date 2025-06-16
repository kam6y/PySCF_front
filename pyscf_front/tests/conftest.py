"""
pytest設定とフィクスチャ
"""
import os
import tempfile
import pytest
from unittest.mock import Mock, patch
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker, Session

# テスト用の環境変数設定
os.environ.update({
    'DB_HOST': 'localhost',
    'DB_PORT': '3306',
    'DB_NAME': 'pyscf_front_test',
    'DB_USER': 'test_user',
    'DB_PASSWORD': 'test_password',
    'PYSCF_FRONT_DEBUG': 'true',
    'PYTEST_QT_API': 'pyside6'
})

import sys
from pathlib import Path
project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))

from pyscf_front.database.models import Base
from pyscf_front.core.molecule import Molecule, Atom


@pytest.fixture
def test_engine():
    """テスト用SQLiteエンジン（各テストで独立）"""
    # 各テストで独立したインメモリーSQLiteデータベースを使用
    engine = create_engine("sqlite:///:memory:", echo=False)
    Base.metadata.create_all(engine)
    return engine


@pytest.fixture
def test_session(test_engine):
    """テスト用セッション（各テストで独立）"""
    SessionLocal = sessionmaker(bind=test_engine)
    session = SessionLocal()
    
    try:
        yield session
    finally:
        session.rollback()
        session.close()
        # セッション終了後にエンジンもクリア
        test_engine.dispose()


@pytest.fixture
def mock_pyscf():
    """PySCFのモック"""
    with patch('pyscf_front.plugins.methods.pyscf_methods.PYSCF_AVAILABLE', True), \
         patch('pyscf_front.plugins.basis_sets.pyscf_basis_sets.PYSCF_AVAILABLE', True):
        # PySCFモジュールのモック
        mock_gto = Mock()
        mock_scf = Mock()
        mock_dft = Mock()
        
        # 分子オブジェクトのモック
        mock_mol = Mock()
        mock_mol.natm = 3
        mock_mol.charge = 0
        mock_mol.spin = 0
        mock_mol.nelectron = 10
        mock_gto.Mole.return_value = mock_mol
        
        # 計算オブジェクトのモック
        mock_mf = Mock()
        mock_mf.kernel.return_value = -76.123456  # 水分子のエネルギー例
        mock_mf.converged = True
        mock_mf.mo_energy = [-1.2, -0.8, -0.6, -0.4, 0.1, 0.3, 0.5]
        mock_mf.dip_moment.return_value = [0.0, 0.0, 1.85]
        mock_mf.mulliken_pop.return_value = (None, [-0.82, 0.41, 0.41])
        mock_mf.niter = 15
        
        mock_scf.RHF.return_value = mock_mf
        mock_scf.UHF.return_value = mock_mf
        mock_dft.RKS.return_value = mock_mf
        mock_dft.UKS.return_value = mock_mf
        
        with patch.multiple(
            'pyscf_front.plugins.methods.pyscf_methods',
            gto=mock_gto,
            scf=mock_scf,
            dft=mock_dft
        ):
            yield {
                'gto': mock_gto,
                'scf': mock_scf,
                'dft': mock_dft,
                'mol': mock_mol,
                'mf': mock_mf
            }


@pytest.fixture
def sample_molecule():
    """テスト用分子（水分子）"""
    molecule = Molecule("H2O", charge=0, multiplicity=1)
    molecule.atoms = [
        Atom("O", 0.0, 0.0, 0.0),
        Atom("H", 0.0, 0.757, 0.587),
        Atom("H", 0.0, -0.757, 0.587)
    ]
    return molecule


@pytest.fixture
def sample_molecule_methane():
    """テスト用分子（メタン）"""
    molecule = Molecule("CH4", charge=0, multiplicity=1)
    molecule.atoms = [
        Atom("C", 0.0, 0.0, 0.0),
        Atom("H", 1.089, 0.0, 0.0),
        Atom("H", -0.363, 1.027, 0.0),
        Atom("H", -0.363, -0.514, 0.890),
        Atom("H", -0.363, -0.514, -0.890)
    ]
    return molecule


@pytest.fixture
def temp_file():
    """一時ファイル"""
    fd, path = tempfile.mkstemp()
    yield path
    os.close(fd)
    os.unlink(path)


@pytest.fixture
def mock_qt_application():
    """QApplicationのモック"""
    with patch('PySide6.QtWidgets.QApplication') as mock_app:
        yield mock_app


@pytest.fixture
def mock_thread_pool():
    """QThreadPoolのモック"""
    with patch('PySide6.QtCore.QThreadPool') as mock_pool:
        mock_instance = Mock()
        mock_pool.return_value = mock_instance
        yield mock_instance