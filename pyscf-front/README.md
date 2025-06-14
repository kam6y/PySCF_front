# PySCF_Front

PySCF_Frontは、量子化学計算を行うためのクロスプラットフォーム対応ネイティブアプリケーションです。研究者や学生が専門的なプログラミング知識なしに、直感的なGUIを通じて高度な量子化学計算を実行できることを目的としています。

## 特徴

- **直感的なGUI**: Flutterを使用したクロスプラットフォーム対応
- **高性能計算**: PySCFによる本格的な量子化学計算
- **データベース管理**: MySQLによる計算結果の体系的管理
- **拡張可能**: プラグインアーキテクチャによる機能拡張
- **MCP連携**: Claude Desktopなどとの連携（オプション）

## 技術スタック

### バックエンド
- Python 3.13
- PySCF (量子化学計算エンジン)
- FastAPI (Web API)
- gRPC (Flutter-Python間通信)
- SQLAlchemy + MySQL (データベース)

### フロントエンド
- Flutter 3.32
- Dart
- Material Design 3

### インフラ
- Docker & Docker Compose
- MySQL 8.0

## セットアップ

### 前提条件

- Python 3.13 (pyenv推奨)
- Flutter 3.32+
- Docker Desktop
- Git

### 1. リポジトリのクローン

```bash
git clone <repository-url>
cd pyscf-front
```

### 2. 環境設定

```bash
# 環境変数ファイルをコピー
cp .env.example .env

# 必要に応じて .env を編集
```

### 3. Docker環境での起動

```bash
# 開発環境を起動
make docker-up-dev

# または直接 docker compose を使用
docker compose -f docker-compose.yml -f docker-compose.dev.yml up -d
```

### 4. ローカル開発環境での起動

```bash
# 開発環境セットアップ
make setup

# バックエンドを起動
make run-backend

# フロントエンドを起動（別ターミナル）
make run-frontend
```

## 開発

### バックエンドの開発

```bash
cd backend

# 仮想環境の作成
python -m venv venv
source venv/bin/activate  # macOS/Linux
# または venv\Scripts\activate  # Windows

# 依存関係のインストール
pip install -r requirements.txt
pip install -e .

# テストの実行
pytest

# APIサーバーの起動
uvicorn src.main:app --reload
```

### フロントエンドの開発

```bash
cd frontend

# 依存関係のインストール
flutter pub get

# アプリケーションの起動
flutter run -d macos  # macOS
flutter run -d windows  # Windows
flutter run -d linux  # Linux

# テストの実行
flutter test
```

## Docker開発環境

### コンテナの起動

```bash
# MySQL + バックエンド + phpMyAdmin
docker compose up -d

# 開発用設定で起動
docker compose -f docker-compose.yml -f docker-compose.dev.yml up -d
```

### アクセス情報

- **バックエンドAPI**: http://localhost:8000
- **phpMyAdmin**: http://localhost:8080
- **MySQL**: localhost:3306
  - ユーザー: `pyscf_user`
  - パスワード: `pyscf_password`
  - データベース: `pyscf_dev`

## MCP連携（オプション）

Claude Desktopとの連携機能を有効にする場合：

1. 設定ファイルでMCPサーバーを有効化
2. Claude Desktopの設定を更新

詳細は [mcp-server-setup.md](mcp-server-setup.md) を参照してください。

## プロジェクト構造

```
pyscf-front/
├── backend/              # Pythonバックエンド
│   ├── src/
│   │   ├── api/         # gRPC/REST API
│   │   ├── core/        # PySCF設定・計算エンジン
│   │   ├── plugins/     # 計算手法プラグイン
│   │   └── utils/       # ユーティリティ
│   ├── tests/
│   └── requirements.txt
├── frontend/             # Flutterフロントエンド
│   ├── lib/
│   │   ├── models/      # データモデル
│   │   ├── views/       # UI画面
│   │   ├── controllers/ # ビジネスロジック
│   │   └── services/    # gRPC通信
│   └── pubspec.yaml
├── database/
│   ├── schema/          # スキーマ定義
│   └── migrations/      # マイグレーション
├── docker/              # Docker設定
├── protos/              # gRPC定義
└── mcp_server/          # MCPサーバー（オプション）
```

## よく使うコマンド

```bash
# 開発環境の完全セットアップ
make setup

# Docker環境の起動/停止
make docker-up
make docker-down

# テストの実行
make test

# キャッシュのクリア
make clean

# ヘルプの表示
make help
```

## ライセンス

Apache License 2.0

## 貢献

プルリクエストやイシューをお待ちしています。詳細は CONTRIBUTING.md をご覧ください。