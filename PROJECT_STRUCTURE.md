# PySCF_native_app プロジェクト構成

## 概要

PySCF_native_appは、PySide6 (Qt 6.7+) を使用してPythonで完全に構築された高度な量子化学計算デスクトップアプリケーションです。PySCFを計算エンジンとして使用し、分子の可視化、計算管理、結果分析を統合したネイティブGUIアプリケーションです。

## アーキテクチャ概要

### 統合Python アーキテクチャ
- **単一言語**: 全コンポーネントがPythonで記述され、シームレスな統合を実現
- **直接統合**: gRPC不要 - GUIと計算エンジン間の直接関数呼び出し
- **ネイティブ性能**: Qt基盤のGUIでOpenGL/VTKアクセラレーション
- **科学計算エコシステム**: NumPy、SciPy、Matplotlib、RDKitとの直接統合

### 技術スタック
- **フロントエンド**: PySide6 6.7.0 (Qt基盤のネイティブGUI)
- **計算エンジン**: PySCF 2.9.0 (量子化学計算)
- **分子処理**: RDKit 2024.03.3 (SMILES処理、3D生成)
- **データベース**: SQLAlchemy 2.0.30 (MySQL/SQLite対応)
- **可視化**: VTK 9.3.0、Matplotlib 3.9.0
- **設定管理**: YAML設定ファイル + 環境変数オーバーライド
- **テスト**: Pytest + Qt GUIテスト対応

## ディレクトリ構成

```
PySCF_native_app/
├── README.md                    # プロジェクト基本説明
├── CLAUDE.md                    # Claude Code用プロジェクト指示書
├── IMPLEMENTATION_STATUS.md     # 実装状況レポート (95%完成)
├── SETUP_COMPLETE.md           # セットアップ完了ドキュメント
├── PROJECT_STRUCTURE.md        # 本ファイル - プロジェクト構成説明
├── mcp-server-setup.md         # MCP サーバーセットアップガイド
├── setup.py                    # Python パッケージ配布設定
├── run_dev.py                  # 開発用起動スクリプト
├── test_environment.py         # 環境検証スクリプト
├── requirements.txt            # 本番依存関係
├── requirements-dev.txt        # 開発依存関係
├── venv/                       # Python仮想環境
│
└── pyscf_front/               # メインアプリケーション
    ├── __init__.py
    ├── main.py                # アプリケーション エントリーポイント
    │
    ├── core/                  # 計算エンジン・コアロジック
    │   ├── __init__.py
    │   ├── molecule.py        # 分子データ管理 (Atom, Molecule, MoleculeBuilder)
    │   ├── calculation_engine.py        # PySCF統合、非同期ジョブ管理
    │   └── calculation_engine_integrated.py  # 統合計算エンジン
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
    │   └── resources/         # UIリソース
    │       ├── icons/         # アプリケーションアイコン
    │       ├── styles/        # スタイルシート
    │       │   └── dark_theme.qss       # ダークテーマ
    │       └── ui_files/      # Qt Designer ファイル
    │
    ├── database/              # データ永続化層
    │   ├── __init__.py
    │   ├── models.py          # SQLAlchemy ORM モデル
    │   ├── repository.py      # データアクセス層 (DAO パターン)
    │   ├── connection.py      # データベース接続管理
    │   └── migrations/        # データベースマイグレーション
    │
    ├── services/              # ビジネスロジック サービス層
    │   ├── __init__.py
    │   ├── calculation_service.py       # 計算管理サービス
    │   ├── instance_service.py          # プロジェクトインスタンス管理
    │   └── molecule_service.py          # 分子操作サービス
    │
    ├── plugins/               # プラグインシステム
    │   ├── __init__.py
    │   ├── methods/           # 計算手法プラグイン
    │   ├── basis_sets/        # 基底関数プラグイン
    │   └── analysis/          # 解析ツールプラグイン
    │
    ├── utils/                 # ユーティリティ関数
    │   ├── __init__.py
    │   ├── config.py          # 設定管理 (YAML + 環境変数)
    │   └── logger.py          # ログ設定・管理
    │
    ├── tests/                 # テストスイート
    │   ├── __init__.py
    │   ├── conftest.py        # Pytest設定・フィクスチャ
    │   ├── test_calculation_engine_integration.py # 統合テスト
    │   └── test_services.py   # サービス層テスト
    │
    ├── translations/          # 国際化対応
    │   └── (翻訳ファイル)
    │
    ├── mcp_server/            # MCP サーバー (オプション)
    │   └── (Claude Desktop統合用)
    │
    └── resources/             # 追加リソース
        └── (追加リソースファイル)
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

#### `calculation_engine.py`
- **CalculationJob データクラス**: ジョブ追跡
- **CalculationEngine**: PySCF 統合計算エンジン
  - スレッド化された計算実行
  - 進捗シグナル対応
  - HF、DFT (B3LYP、PBE)、post-HF 手法対応
  - 非同期ジョブ管理

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

### 3. データ管理層 (`database/`)

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
- データベース接続管理
- MySQL/SQLite 対応

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

### 5. 設定管理 (`utils/`)

#### `config.py`
- YAML 設定ファイル管理
- 環境変数オーバーライド対応
- **設定セクション**:
  - データベース設定
  - アプリケーション設定
  - 計算エンジン設定
  - GUI 設定
  - MCP サーバー設定

### 6. テスト (`tests/`)

#### `conftest.py`
- Pytest 設定・フィクスチャ
- テスト用データベース (SQLite インメモリ)
- PySCF モック
- サンプル分子 (水、メタン)
- Qt アプリケーションテスト対応

#### テストクラス
- **統合テスト**: 計算エンジンとデータベース統合
- **サービス層テスト**: ビジネスロジックテスト
- **GUI テスト**: pytest-qt 使用

## データフロー

1. **ユーザー入力** → PySide6 GUI ダイアログ
2. **GUI** → **サービス層** (直接 Python 関数呼び出し)
3. **サービス層** → **コア計算エンジン** (PySCF)
4. **計算実行** → QThreadPool でバックグラウンド処理
5. **結果** → SQLAlchemy 経由でデータベース保存
6. **可視化** → VTK/OpenGL で 3D 分子構造レンダリング

## 実装状況 (95% 完成)

### ✅ 完了済み
- 完全な分子入力系統 (SMILES、XYZ、プリセット)
- 完全な PySCF 計算エンジン統合
- 包括的な GUI とジョブ管理
- 非同期計算実行と進捗追跡
- 基本的な 3D 分子表示
- データベース永続化とサービス層統合

### 🔄 進行中
- 高度な 3D 可視化 (60% 完成)
- 追加の解析ツール

## 外部統合

### MCP サーバー (オプション)
- Model Context Protocol サーバー
- Claude Desktop との統合
- 自然言語による計算操作
- セキュリティ制限 (IP制限、レート制限)

## 開発環境

### 依存関係管理
- **requirements.txt**: 本番環境依存関係
- **requirements-dev.txt**: 開発環境依存関係
- **venv/**: Python 仮想環境

### 開発ツール
- **Black**: コードフォーマッター (88文字、Python 3.13対応)
- **MyPy**: 型チェック
- **Flake8**: PEP 8 コンプライアンス
- **Pytest**: テストフレームワーク

### IDE サポート
- VS Code 設定完備
- Python インタープリター自動設定
- Qt ツール統合 (Designer、UIC、RCC)
- デバッグ・実行設定

## 使用方法

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
pytest pyscf_front/tests/ -v
```

このプロジェクトは、量子化学計算のための成熟した、よく設計されたデスクトップアプリケーションであり、プロフェッショナルな開発環境と包括的なテスト、設定管理、外部統合機能を提供します。