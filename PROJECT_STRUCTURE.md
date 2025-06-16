# PySCF_native_app プロジェクト構成

## 概要

PySCF_native_appは、PySide6 (Qt 6.7+) を使用してPythonで完全に構築された高度な量子化学計算デスクトップアプリケーションです。PySCFを計算エンジンとして使用し、分子の可視化、計算管理、結果分析を統合したネイティブGUIアプリケーション、そして拡張可能なプラグインシステムを提供します。

## アーキテクチャ概要

### 統合Python アーキテクチャ
- **単一言語**: 全コンポーネントがPythonで記述され、シームレスな統合を実現
- **直接統合**: gRPC不要 - GUIと計算エンジン間の直接関数呼び出し
- **ネイティブ性能**: Qt基盤のGUIでOpenGL/VTKアクセラレーション
- **科学計算エコシステム**: NumPy、SciPy、Matplotlib、RDKitとの直接統合
- **プラグインアーキテクチャ**: 拡張可能な計算手法、基底関数、解析ツールのプラグインシステム

### 技術スタック
- **フロントエンド**: PySide6 6.7.0+ (Qt基盤のネイティブGUI)
- **計算エンジン**: PySCF 2.9.0+ (量子化学計算) + 統合プラグインシステム
- **分子処理**: RDKit 2024.03.3+ (SMILES処理、3D生成)
- **データベース**: SQLAlchemy 2.0.30+ (SQLite専用)
- **可視化**: VTK 9.3.0+、Matplotlib 3.9.0+
- **設定管理**: 環境変数ベース設定システム
- **テスト**: Pytest + pytest-qt (GUI テスト対応)

## ディレクトリ構成

```
PySCF_native_app/
├── README.md                    # プロジェクト基本説明
├── CLAUDE.md                    # Claude Code用プロジェクト指示書
├── IMPLEMENTATION_STATUS.md     # 実装状況レポート
├── SETUP_COMPLETE.md           # セットアップ完了ドキュメント
├── PROJECT_STRUCTURE.md        # 本ファイル - プロジェクト構成説明
├── LICENSE                     # ライセンス情報
├── mcp-server-setup.md         # MCP サーバーセットアップガイド
├── setup.py                    # Python パッケージ配布設定
├── run_dev.py                  # 開発用起動スクリプト
├── test_environment.py         # 環境検証スクリプト
├── requirements.txt            # 本番依存関係
├── requirements-dev.txt        # 開発依存関係
├── .env                        # 環境変数設定ファイル
├── .env.example               # 環境変数設定例
├── .gitignore                 # Git除外設定
├── venv/                      # Python仮想環境
├── data/                      # データディレクトリ
│   └── pyscf_front.db         # SQLiteデータベースファイル
├── logs/                      # ログディレクトリ
│
└── pyscf_front/               # メインアプリケーション
    ├── __init__.py
    ├── main.py                # アプリケーション エントリーポイント
    │
    ├── core/                  # 計算エンジン・コアロジック
    │   ├── __init__.py
    │   ├── molecule.py        # 分子データ管理 (Atom, Molecule, MoleculeBuilder)
    │   ├── calculation_engine.py         # レガシー計算エンジン
    │   └── calculation_engine_unified.py # 統合計算エンジン（プラグイン対応）
    │
    ├── gui/                   # PySide6 GUI コンポーネント
    │   ├── __init__.py
    │   ├── main_window.py     # メインウィンドウ (ドック可能パネル)
    │   ├── dialogs/           # モーダルダイアログ
    │   │   ├── __init__.py
    │   │   ├── calculation_dialog.py    # 計算設定ダイアログ
    │   │   └── molecule_input_dialog.py # 分子入力ダイアログ
    │   ├── widgets/           # 再利用可能UIウィジェット
    │   │   ├── __init__.py
    │   │   └── molecule_viewer.py       # 分子管理ウィジェット
    │   └── resources/         # UIリソース (空のディレクトリ)
    │
    ├── database/              # データ永続化層（SQLite専用）
    │   ├── __init__.py
    │   ├── models.py          # SQLAlchemy ORM モデル
    │   ├── repository.py      # データアクセス層 (DAO パターン)
    │   ├── connection.py      # SQLiteデータベース接続管理
    │   └── migrations/        # データベースマイグレーション (空)
    │
    ├── services/              # ビジネスロジック サービス層
    │   ├── __init__.py
    │   ├── calculation_service.py       # 計算管理サービス
    │   ├── instance_service.py          # プロジェクトインスタンス管理
    │   └── molecule_service.py          # 分子操作サービス
    │
    ├── plugins/               # プラグインシステム
    │   ├── __init__.py
    │   ├── base.py            # プラグイン基底クラス・管理システム
    │   ├── methods/           # 計算手法プラグイン
    │   │   ├── __init__.py
    │   │   └── pyscf_methods.py        # PySCF手法プラグイン
    │   ├── basis_sets/        # 基底関数プラグイン
    │   │   ├── __init__.py
    │   │   └── pyscf_basis_sets.py     # PySCF基底関数プラグイン
    │   └── analysis/          # 解析ツールプラグイン (空)
    │
    ├── utils/                 # ユーティリティ関数
    │   ├── __init__.py
    │   ├── config.py          # 設定管理 (環境変数ベース)
    │   └── logger.py          # ログ設定・管理
    │
    ├── tests/                 # テストスイート
    │   ├── __init__.py
    │   ├── conftest.py        # Pytest設定・フィクスチャ
    │   ├── test_gui_components.py       # GUI コンポーネントテスト
    │   ├── test_plugin_system.py       # プラグインシステムテスト
    │   ├── test_services.py            # サービス層テスト
    │   └── test_unified_engine.py      # 統合エンジンテスト
    │
    ├── translations/          # 国際化対応 (空)
    │   └── (翻訳ファイル)
    │
    ├── mcp_server/            # MCP サーバー (オプション、空)
    │   └── (Claude Desktop統合用)
    │
    ├── docs/                  # ドキュメント (空)
    │   └── (プロジェクト文書)
    │
    └── resources/             # アプリケーションリソース
        ├── icons/             # アプリケーションアイコン (空)
        ├── styles/            # スタイルシート
        │   └── dark_theme.qss # ダークテーマ
        └── ui_files/          # Qt Designer ファイル (空)
```

## 主要コンポーネント詳細

### 1. コア計算エンジン (`core/`)

#### `molecule.py`
- **Atom データクラス**: 原子データ (記号、座標、電荷)
- **Molecule クラス**: 分子データ管理
  - XYZ ファイル解析
  - 分子式自動生成
  - 重心計算
  - RDKit 統合による 3D 構造生成
- **MoleculeBuilder**: 分子構築ヘルパー
  - SMILES 文字列からの分子生成
  - 水、メタン、ベンゼンなどのプリセット分子
- **MoleculeManager**: 分子コレクション管理

#### `calculation_engine.py` (レガシー)
- **CalculationJob データクラス**: ジョブ追跡
- **CalculationEngine**: 従来のPySCF 統合計算エンジン
  - スレッド化された計算実行
  - 進捗シグナル対応
  - HF、DFT (B3LYP、PBE)、post-HF 手法対応

#### `calculation_engine_unified.py` (現在の主要エンジン)
- **UnifiedCalculationEngine**: プラグインシステム対応の統合計算エンジン
  - SQLite専用データベース永続化
  - プラグインベースの計算手法・基底関数管理
  - 非同期ジョブ管理と進捗追跡
  - 自動的な計算時間推定
  - 結果の詳細解析（HOMO/LUMO、双極子モーメント、Mulliken電荷など）

### 2. GUI コンポーネント (`gui/`)

#### `main_window.py`
- メインアプリケーションウィンドウ
- ドック可能パネル (分子管理、計算設定、結果表示)
- 分子マネージャーとの統合
- 計算エンジン協調
- メニューバー・ステータスバー

#### `dialogs/`
- **calculation_dialog.py**: 計算設定ダイアログ
- **molecule_input_dialog.py**: 分子入力インターフェース

#### `widgets/`
- **molecule_viewer.py**: 分子管理ウィジェット
- 再利用可能な GUI コンポーネント

### 3. データ管理層 (`database/`) - SQLite専用

#### `models.py`
- SQLAlchemy ORM モデル定義
- **Enum クラス**: InstanceStatus、CalculationStatus、JobStatus
- **データベーススキーマ**:
  - `instances`: プロジェクトコンテナ
  - `molecules`: 分子構造データ
  - `calculations`: 計算ジョブ情報
  - `results`: 計算結果
  - `job_queue`: ジョブキュー管理

#### `repository.py`
- データアクセス層 (DAO パターン)
- **リポジトリクラス**:
  - InstanceRepository: インスタンス操作
  - MoleculeRepository: 分子データ操作
  - CalculationRepository: 計算管理
  - ResultRepository: 結果管理
  - JobQueueRepository: ジョブキュー操作

#### `connection.py`
- **DatabaseConfig**: SQLite専用設定クラス
- **DatabaseManager**: SQLiteデータベース接続管理
- ローカルファイルベースのデータベース (`data/pyscf_front.db`)
- 自動テーブル作成とスキーマ管理

### 4. サービス層 (`services/`)

#### `calculation_service.py`
- 計算ビジネスロジック
- 計算ライフサイクル管理
- 結果保存とキュー管理

#### `instance_service.py`
- プロジェクトインスタンス管理
- 完全なインスタンス作成 (分子 + 計算設定)

#### `molecule_service.py`
- 分子操作サービス
- データベース統合
- SMILES/XYZ からの分子作成

### 5. プラグインシステム (`plugins/`)

#### `base.py`
- **PluginInterface**: 全プラグインの基底クラス
- **CalculationMethodPlugin**: 計算手法プラグインの基底クラス
- **BasisSetPlugin**: 基底関数プラグインの基底クラス
- **AnalysisPlugin**: 解析プラグインの基底クラス
- **PluginManager**: プラグインの登録・管理・実行制御

#### `methods/pyscf_methods.py`
- **PySCFMethodPlugin**: PySCF対応計算手法プラグイン
- サポート手法: HF、B3LYP、PBE、M06、wB97X-D など
- 計算時間推定とパラメータ検証機能

#### `basis_sets/pyscf_basis_sets.py`
- **PySCFBasisSetPlugin**: PySCF対応基底関数プラグイン
- サポート基底関数: STO-3G、6-31G系、cc-pVDZ系、def2系 など
- 原子種との互換性チェック機能

### 6. 設定管理 (`utils/`)

#### `config.py`
- 環境変数ベース設定管理
- **主要設定項目**:
  - データベース設定 (`DB_PATH`)
  - アプリケーション設定 (`PYSCF_FRONT_DEBUG`)
  - ログレベル設定 (`PYSCF_FRONT_LOG_LEVEL`)
  - GUI設定 (`PYSCF_FRONT_THEME`, `PYSCF_FRONT_LANGUAGE`)

### 7. テスト (`tests/`)

#### `conftest.py`
- Pytest 設定・フィクスチャ
- テスト用データベース (SQLite インメモリ)
- PySCF モック
- サンプル分子 (水、メタン)
- Qt アプリケーションテスト対応

#### テストモジュール
- **test_gui_components.py**: GUI コンポーネントテスト（pytest-qt使用）
- **test_plugin_system.py**: プラグインシステムの登録・実行テスト
- **test_services.py**: サービス層のビジネスロジックテスト
- **test_unified_engine.py**: 統合計算エンジンの機能テスト

## データフロー（現在のアーキテクチャ）

1. **ユーザー入力** → PySide6 GUI ダイアログ・ウィジェット
2. **GUI層** → **サービス層** (直接 Python 関数呼び出し)
3. **サービス層** → **統合計算エンジン** (プラグインシステム経由)
4. **プラグイン選択** → 手法・基底関数プラグインの動的ロード
5. **計算実行** → QThreadPool でバックグラウンド処理 + PySCF計算
6. **進捗更新** → Qt シグナル/スロットによる GUI 更新
7. **結果処理** → SQLAlchemy 経由で SQLite データベース保存
8. **可視化** → VTK/OpenGL で 3D 分子構造レンダリング（予定）

## 実装状況

### ✅ 完了済み
- **コアシステム**: 分子入力・管理システム (SMILES、XYZ、プリセット)
- **計算エンジン**: 統合計算エンジンとプラグインシステム
- **データベース**: SQLite専用データ永続化システム
- **GUI基盤**: PySide6ベースのメインウィンドウとダイアログ
- **サービス層**: 完全なビジネスロジック実装
- **テスト**: 包括的なユニット・統合テストスイート

### 🔄 進行中・予定
- **3D可視化**: VTK統合による高度な分子表示
- **解析ツール**: 追加の結果解析プラグイン
- **MCP サーバー**: Claude Desktop統合機能
- **国際化**: 多言語対応システム

## 外部統合

### MCP サーバー (オプション)
- Model Context Protocol サーバー
- Claude Desktop との統合
- 自然言語による計算操作
- セキュリティ制限 (IP制限、レート制限)

## 開発環境

### 依存関係管理
- **requirements.txt**: 本番環境依存関係 (PySCF、PySide6、SQLAlchemy等)
- **requirements-dev.txt**: 開発環境依存関係 (pytest、Black、MyPy等)
- **venv/**: Python 3.13 仮想環境

### 開発ツール設定
- **Black**: コードフォーマッター (88文字、Python 3.13対応)
- **MyPy**: 静的型チェック
- **Flake8**: PEP 8 コンプライアンス
- **Pytest**: テストフレームワーク + pytest-qt (GUI テスト)

### 環境設定
- **.env**: 環境変数設定ファイル (SQLiteパスなど)
- **.env.example**: 環境変数設定テンプレート
- **data/**: SQLiteデータベースファイル格納ディレクトリ
- **logs/**: アプリケーションログ格納ディレクトリ

### IDE サポート
- VS Code 設定完備 (推奨)
- Python インタープリター自動設定 (venv/bin/python)
- Qt ツール統合 (Designer、UIC、RCC)
- デバッグ・実行設定

## 使用方法

### 初期セットアップ
```bash
# 仮想環境作成・有効化
python -m venv venv
source venv/bin/activate  # Windows: venv\Scripts\activate

# 依存関係インストール
pip install -r requirements.txt
pip install -r requirements-dev.txt

# 環境設定
cp .env.example .env  # 必要に応じて編集
```

### 開発モード起動
```bash
python run_dev.py
```

### 環境設定確認
```bash
python test_environment.py
```

### テスト実行
```bash
# 全テスト実行
pytest pyscf_front/tests/ -v

# GUI テスト実行
pytest pyscf_front/tests/test_gui_components.py --qt-test-timeout=5000

# カバレッジ付きテスト
pytest pyscf_front/tests/ --cov=pyscf_front
```

### コード品質チェック
```bash
# フォーマット
black pyscf_front/ --line-length=88

# 型チェック
mypy pyscf_front/

# リント
flake8 pyscf_front/
```

## プロジェクトの特徴

PySCF_native_appは、以下の特徴を持つ成熟した量子化学計算デスクトップアプリケーションです：

- **統合アーキテクチャ**: 単一Python言語による seamless な統合
- **モジュラー設計**: プラグインシステムによる高い拡張性
- **データ永続化**: SQLite による軽量で高速なローカルデータベース
- **プロフェッショナル品質**: 包括的なテスト、型注釈、コード品質管理
- **ユーザーフレンドリー**: 直感的なQt/PySide6 GUI
- **科学計算特化**: PySCF統合による高精度量子化学計算