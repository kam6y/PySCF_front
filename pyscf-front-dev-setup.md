# PySCF_Front 開発環境構築ガイド v1.0

## 目次

1. [概要](#1-概要)
2. [前提条件](#2-前提条件)
3. [開発環境の全体構成](#3-開発環境の全体構成)
4. [Python環境構築](#4-python環境構築)
5. [Flutter環境構築](#5-flutter環境構築)
6. [データベース環境構築](#6-データベース環境構築)
7. [バックエンド開発環境構築](#7-バックエンド開発環境構築)
8. [フロントエンド開発環境構築](#8-フロントエンド開発環境構築)
9. [GPU開発環境構築（オプション）](#9-gpu開発環境構築オプション)
10. [MCPサーバー開発環境（オプション）](#10-mcpサーバー開発環境オプション)
11. [統合開発環境セットアップ](#11-統合開発環境セットアップ)
12. [トラブルシューティング](#12-トラブルシューティング)

---

## 1. 概要

このガイドは、PySCF_Frontの開発環境を構築するための詳細な手順を提供します。要件定義書に基づき、最新の技術スタックを使用した開発環境を構築します。

### 1.1 注意事項

**重要**: Python 3.14は2025年6月現在ベータ版（3.14.0b2）であり、一部のライブラリが未対応の可能性があります。本ガイドでは以下の2つのアプローチを提供します：

1. **推奨**: Python 3.13.0を使用した安定版環境
2. **実験的**: Python 3.14.0b2を使用した最新版環境

## 2. 前提条件

### 2.1 ハードウェア要件

**最小要件**
- CPU: 4コア以上（Intel i5 または AMD Ryzen 5 相当）
- RAM: 16GB以上
- ストレージ: 50GB以上の空き容量（SSD推奨）
- GPU: なし（CPU開発のみ）

**推奨要件**
- CPU: 8コア以上（Intel i7 または AMD Ryzen 7 相当）
- RAM: 32GB以上
- ストレージ: 200GB以上の空き容量（NVMe SSD推奨）
- GPU: NVIDIA RTX 3060以上（CUDA 11.0+対応）

### 2.2 OS要件

| OS | バージョン | 優先度 |
|----|-----------|--------|
| macOS | 13.0 (Ventura) 以上 | 高（開発優先） |
| Windows | 10 (64bit) / 11 | 高（GPU開発） |
| Ubuntu | 22.04 LTS / 24.04 LTS | 中 |

### 2.3 必須ソフトウェア

- Git 2.40+
- Docker Desktop 24.0+
- Visual Studio Code（推奨）またはその他のIDE

## 3. 開発環境の全体構成

```
pyscf-front/
├── backend/            # Pythonバックエンド
│   ├── src/
│   ├── tests/
│   ├── requirements.txt
│   └── pyproject.toml
├── frontend/           # Flutterフロントエンド
│   ├── lib/
│   ├── test/
│   ├── pubspec.yaml
│   └── analysis_options.yaml
├── database/           # データベース設定
│   ├── schema/
│   └── migrations/
├── docker/             # Docker設定
│   ├── backend/
│   ├── frontend/
│   └── database/
├── mcp_server/         # MCPサーバー（オプション）
│   ├── src/
│   └── requirements.txt
├── protos/             # gRPC定義
│   └── calculation.proto
├── docker-compose.yml
├── docker-compose.dev.yml
└── README.md
```

## 4. Python環境構築

### 4.1 pyenvのインストール

**macOS**
```bash
# Homebrewを使用
brew install pyenv
echo 'export PYENV_ROOT="$HOME/.pyenv"' >> ~/.zshrc
echo 'command -v pyenv >/dev/null || export PATH="$PYENV_ROOT/bin:$PATH"' >> ~/.zshrc
echo 'eval "$(pyenv init -)"' >> ~/.zshrc
source ~/.zshrc
```

**Ubuntu/WSL**
```bash
# 依存パッケージのインストール
sudo apt update
sudo apt install -y make build-essential libssl-dev zlib1g-dev \
libbz2-dev libreadline-dev libsqlite3-dev wget curl llvm libncurses5-dev \
libncursesw5-dev xz-utils tk-dev libffi-dev liblzma-dev python3-openssl git

# pyenvのインストール
curl https://pyenv.run | bash

# .bashrcに追加
echo 'export PYENV_ROOT="$HOME/.pyenv"' >> ~/.bashrc
echo 'command -v pyenv >/dev/null || export PATH="$PYENV_ROOT/bin:$PATH"' >> ~/.bashrc
echo 'eval "$(pyenv init -)"' >> ~/.bashrc
source ~/.bashrc
```

**Windows**
PowerShellを管理者権限で実行：
```powershell
# pyenv-winのインストール
git clone https://github.com/pyenv-win/pyenv-win.git "$HOME\.pyenv"

# 環境変数の設定
[System.Environment]::SetEnvironmentVariable('PYENV', "$HOME\.pyenv\pyenv-win", "User")
[System.Environment]::SetEnvironmentVariable('PYENV_ROOT', "$HOME\.pyenv\pyenv-win", "User")
[System.Environment]::SetEnvironmentVariable('PYENV_HOME', "$HOME\.pyenv\pyenv-win", "User")

# PATHに追加
$path = [System.Environment]::GetEnvironmentVariable('PATH', "User")
[System.Environment]::SetEnvironmentVariable('PATH', "$HOME\.pyenv\pyenv-win\bin;$HOME\.pyenv\pyenv-win\shims;$path", "User")
```

### 4.2 Pythonのインストール

#### オプション1: 安定版環境（推奨）
```bash
# Python 3.13.0のインストール
pyenv install 3.13.0
pyenv global 3.13.0

# 確認
python --version  # Python 3.13.0
```

#### オプション2: 実験的環境（Python 3.14）
```bash
# Python 3.14.0b2のインストール
pyenv install 3.14.0b2
pyenv global 3.14.0b2

# 確認
python --version  # Python 3.14.0b2
```

### 4.3 仮想環境の作成

```bash
# プロジェクトディレクトリの作成
mkdir -p ~/projects/pyscf-front
cd ~/projects/pyscf-front

# 仮想環境の作成
python -m venv venv

# 仮想環境の有効化
# macOS/Linux
source venv/bin/activate

# Windows
venv\Scripts\activate

# pipのアップグレード
pip install --upgrade pip
```

## 5. Flutter環境構築

### 5.1 Flutter SDKのインストール

**macOS**
```bash
# Homebrewを使用（推奨）
brew install flutter

# または手動インストール
git clone https://github.com/flutter/flutter.git -b stable ~/flutter
echo 'export PATH="$PATH:$HOME/flutter/bin"' >> ~/.zshrc
source ~/.zshrc
```

**Ubuntu/WSL**
```bash
# 必要なパッケージのインストール
sudo apt-get install -y curl git unzip xz-utils zip libglu1-mesa

# Flutterのダウンロード
cd ~
git clone https://github.com/flutter/flutter.git -b stable

# PATHの設定
echo 'export PATH="$PATH:$HOME/flutter/bin"' >> ~/.bashrc
source ~/.bashrc
```

**Windows**
1. [Flutter公式サイト](https://docs.flutter.dev/get-started/install/windows)から最新版をダウンロード
2. C:\flutter に展開
3. システム環境変数のPATHに `C:\flutter\bin` を追加

### 5.2 Flutter環境の確認

```bash
# Flutter Doctorの実行
flutter doctor -v

# 必要に応じて追加設定
flutter config --enable-windows-desktop  # Windows
flutter config --enable-macos-desktop    # macOS
flutter config --enable-linux-desktop    # Linux

# Flutter 3.32へのアップグレード（必要な場合）
flutter upgrade
flutter --version  # Flutter 3.32.0 • channel stable
```

### 5.3 開発用エディタの設定

**VS Code**
```bash
# Flutter/Dart拡張機能のインストール
code --install-extension Dart-Code.dart-code
code --install-extension Dart-Code.flutter
```

**Android Studio**
- Flutter プラグインをインストール
- Dart プラグインをインストール

## 6. データベース環境構築

### 6.1 MySQL 8.0のインストール

**Docker を使用（推奨）**
```yaml
# docker-compose.yml
version: '3.8'
services:
  mysql:
    image: mysql:8.0
    restart: always
    environment:
      MYSQL_ROOT_PASSWORD: dev_root_password
      MYSQL_DATABASE: pyscf_dev
      MYSQL_USER: pyscf_user
      MYSQL_PASSWORD: pyscf_password
    ports:
      - "3306:3306"
    volumes:
      - mysql_data:/var/lib/mysql
      - ./database/schema:/docker-entrypoint-initdb.d
    command: 
      - --character-set-server=utf8mb4
      - --collation-server=utf8mb4_unicode_ci
      - --default-authentication-plugin=mysql_native_password

volumes:
  mysql_data:
```

**ネイティブインストール**

macOS:
```bash
brew install mysql@8.0
brew services start mysql@8.0
```

Ubuntu:
```bash
sudo apt update
sudo apt install mysql-server-8.0
sudo systemctl start mysql
```

Windows:
[MySQL公式サイト](https://dev.mysql.com/downloads/mysql/)からインストーラーをダウンロード

### 6.2 データベースの初期設定

```bash
# MySQLにログイン
mysql -u root -p

# データベースとユーザーの作成
CREATE DATABASE pyscf_dev CHARACTER SET utf8mb4 COLLATE utf8mb4_unicode_ci;
CREATE USER 'pyscf_user'@'localhost' IDENTIFIED BY 'pyscf_password';
GRANT ALL PRIVILEGES ON pyscf_dev.* TO 'pyscf_user'@'localhost';
FLUSH PRIVILEGES;
```

## 7. バックエンド開発環境構築

### 7.1 プロジェクト構造の作成

```bash
cd ~/projects/pyscf-front

# ディレクトリ構造の作成
mkdir -p backend/{src,tests,docs}
mkdir -p backend/src/{api,core,plugins,utils}
mkdir -p backend/src/plugins/{methods,basis_sets,analysis}
```

### 7.2 依存関係のインストール

**requirements.txt の作成**
```txt
# Core dependencies
pyscf==2.9.0
numpy==1.26.4
scipy==1.13.1
h5py==3.11.0

# Database
mysql-connector-python==9.3.0
sqlalchemy==2.0.30
alembic==1.13.1

# gRPC
grpcio==1.73.0
grpcio-tools==1.73.0
protobuf==5.27.0

# Web framework (for API)
fastapi==0.111.0
uvicorn[standard]==0.30.1

# Scientific computing
rdkit==2024.03.3
pandas==2.2.2
matplotlib==3.9.0

# Testing
pytest==8.2.2
pytest-asyncio==0.23.7
pytest-cov==5.0.0

# Code quality
black==24.4.2
isort==5.13.2
flake8==7.0.0
mypy==1.10.0

# Documentation
sphinx==7.3.7
sphinx-rtd-theme==2.0.0

# GPU support (optional)
# gpu4pyscf-cuda11x==1.4.0  # CUDA 11.x
# gpu4pyscf-cuda12x==1.4.0  # CUDA 12.x
```

**Python 3.14用の調整（実験的）**
```txt
# 一部のパッケージは開発版を使用
--pre  # プレリリース版を許可

# 代替パッケージ
# mysql-connector-python の代わりに
pymysql==1.1.1
cryptography==42.0.8

# grpcio のビルドが必要な場合
# ソースからビルドする設定を追加
```

**インストール**
```bash
# 仮想環境が有効化されていることを確認
pip install -r backend/requirements.txt

# 開発用依存関係
pip install -e .
```

### 7.3 PySCF設定

**backend/src/core/pyscf_config.py**
```python
"""PySCF設定モジュール"""
import os
import pyscf
from pyscf import lib

# メモリ設定
lib.num_threads(int(os.environ.get('OMP_NUM_THREADS', 4)))
pyscf.lib.misc.TMPDIR = os.environ.get('PYSCF_TMPDIR', '/tmp')

# 計算精度設定
CONV_TOL = 1e-9
CONV_TOL_GRAD = 1e-6

# GPU設定（オプション）
USE_GPU = os.environ.get('PYSCF_USE_GPU', 'false').lower() == 'true'

if USE_GPU:
    try:
        import gpu4pyscf
        GPU_AVAILABLE = True
    except ImportError:
        GPU_AVAILABLE = False
        print("Warning: GPU4PySCF not available")
else:
    GPU_AVAILABLE = False
```

### 7.4 gRPCサービスの実装

**プロトコル定義のコンパイル**
```bash
# protosディレクトリの作成
mkdir -p protos
cd protos

# calculation.proto の作成（要件定義書から）
cat > calculation.proto << 'EOF'
syntax = "proto3";

package pyscf_front;

// ... (要件定義書のプロトコル定義をここに挿入)
EOF

# Pythonコードの生成
cd ..
python -m grpc_tools.protoc -I./protos --python_out=./backend/src/api --grpc_python_out=./backend/src/api ./protos/calculation.proto
```

### 7.5 プラグインアーキテクチャの実装

**backend/src/plugins/base.py**
```python
"""プラグインベースクラス"""
from abc import ABC, abstractmethod
from typing import Dict, Any, Optional
import numpy as np

class CalculationPlugin(ABC):
    """計算プラグインの基底クラス"""
    
    @classmethod
    @abstractmethod
    def get_name(cls) -> str:
        """プラグイン名を返す"""
        pass
    
    @classmethod
    @abstractmethod
    def get_version(cls) -> str:
        """プラグインバージョンを返す"""
        pass
    
    @abstractmethod
    def validate_input(self, molecule: Any, parameters: Dict[str, Any]) -> bool:
        """入力パラメータの検証"""
        pass
    
    @abstractmethod
    def run_calculation(self, molecule: Any, parameters: Dict[str, Any]) -> Dict[str, Any]:
        """計算の実行"""
        pass
    
    @abstractmethod
    def parse_results(self, raw_output: Any) -> Dict[str, Any]:
        """結果の解析"""
        pass
```

## 8. フロントエンド開発環境構築

### 8.1 Flutterプロジェクトの作成

```bash
cd ~/projects/pyscf-front

# Flutterプロジェクトの作成
flutter create --org com.pyscf --project-name pyscf_front frontend
cd frontend

# 必要なプラットフォームを有効化
flutter config --enable-windows-desktop
flutter config --enable-macos-desktop
flutter config --enable-linux-desktop
```

### 8.2 依存関係の設定

**pubspec.yaml**
```yaml
name: pyscf_front
description: Quantum chemistry calculation native application
version: 1.0.0+1

environment:
  sdk: '>=3.2.0 <4.0.0'
  flutter: '>=3.32.0'

dependencies:
  flutter:
    sdk: flutter
  
  # UI Components
  cupertino_icons: ^1.0.8
  material_design_icons_flutter: ^7.0.7296
  
  # State Management
  provider: ^6.1.2
  riverpod: ^2.5.1
  
  # Navigation
  go_router: ^14.1.4
  
  # gRPC
  grpc: ^4.0.0
  protobuf: ^3.1.0
  
  # 3D Visualization
  flutter_gl: ^0.0.21
  vector_math: ^2.1.4
  
  # File Operations
  file_picker: ^8.0.3
  path_provider: ^2.1.3
  
  # Database
  sqflite: ^2.3.3
  
  # Utilities
  uuid: ^4.4.0
  intl: ^0.19.0
  logger: ^2.3.0
  
  # Charts
  fl_chart: ^0.68.0
  
  # Platform specific
  window_manager: ^0.3.9  # Desktop window management

dev_dependencies:
  flutter_test:
    sdk: flutter
  flutter_lints: ^4.0.0
  build_runner: ^2.4.9
  mockito: ^5.4.4
  integration_test:
    sdk: flutter

flutter:
  uses-material-design: true
  
  assets:
    - assets/images/
    - assets/icons/
    - assets/data/
  
  fonts:
    - family: NotoSansJP
      fonts:
        - asset: assets/fonts/NotoSansJP-Regular.ttf
        - asset: assets/fonts/NotoSansJP-Bold.ttf
          weight: 700
```

### 8.3 プロトコルバッファの設定

```bash
# Flutter用のprotobufコード生成
cd frontend

# protoc-pluginのインストール
dart pub global activate protoc_plugin

# プロトコルファイルからDartコードを生成
protoc --dart_out=lib/generated -I../protos ../protos/calculation.proto
protoc --dart_out=grpc:lib/generated -I../protos ../protos/calculation.proto
```

### 8.4 プロジェクト構造

```bash
# lib/ディレクトリ構造の作成
mkdir -p lib/{models,views,controllers,services,widgets,utils}
mkdir -p lib/views/{home,molecule,calculation,results}
mkdir -p lib/widgets/{common,molecule,calculation}
```

## 9. GPU開発環境構築（オプション）

### 9.1 CUDA環境のセットアップ

**Windows/Linux**
1. [NVIDIA Driver](https://www.nvidia.com/Download/index.aspx)のインストール（最新版）
2. [CUDA Toolkit 11.8](https://developer.nvidia.com/cuda-11-8-0-download-archive)または[CUDA Toolkit 12.x](https://developer.nvidia.com/cuda-downloads)のインストール

```bash
# CUDAバージョンの確認
nvcc --version
nvidia-smi
```

### 9.2 GPU4PySCFのインストール

```bash
# CUDA 11.x用
pip install gpu4pyscf-cuda11x==1.4.0

# CUDA 12.x用
pip install gpu4pyscf-cuda12x==1.4.0

# 開発版（ソースから）
git clone https://github.com/pyscf/gpu4pyscf.git
cd gpu4pyscf
pip install -e .
```

### 9.3 GPU環境のテスト

```python
# test_gpu.py
import pyscf
from pyscf import gto, dft

# GPU利用可能性の確認
try:
    mol = gto.M(
        atom='H 0 0 0; H 0 0 0.74',
        basis='def2-tzvpp'
    )
    mf = dft.RKS(mol, xc='LDA').density_fit()
    
    # GPUへの変換
    mf_gpu = mf.to_gpu()
    e_gpu = mf_gpu.kernel()
    
    print(f"GPU計算成功: Energy = {e_gpu}")
except Exception as e:
    print(f"GPU計算失敗: {e}")
```

## 10. MCPサーバー開発環境（オプション）

### 10.1 MCPサーバーの設定

**mcp_server/requirements.txt**
```txt
# MCP dependencies
anthropic-model-context-protocol==0.1.0
pydantic==2.7.4
fastapi==0.111.0
uvicorn[standard]==0.30.1
websockets==12.0
```

### 10.2 MCPサーバーの実装

**mcp_server/src/server.py**
```python
"""MCPサーバー実装"""
import asyncio
from typing import List, Dict, Any
from mcp import Server, Tool, Resource
from mcp.types import TextContent
import sys
sys.path.append('../../backend/src')

class PySCFMCPServer:
    def __init__(self):
        self.server = Server("pyscf-front")
        self.setup_tools()
    
    def setup_tools(self):
        @self.server.tool("create_molecule")
        async def create_molecule(name: str, description: str) -> Dict[str, Any]:
            """分子構造を生成"""
            # 実装
            return {"status": "success", "molecule_id": "mol_123"}
        
        @self.server.tool("run_calculation")
        async def run_calculation(
            molecule_id: str, 
            method: str = "HF",
            basis: str = "6-31G"
        ) -> Dict[str, Any]:
            """計算を実行"""
            # 実装
            return {"status": "submitted", "job_id": "job_456"}
    
    async def run(self):
        """サーバーの起動"""
        await self.server.run()

if __name__ == "__main__":
    server = PySCFMCPServer()
    asyncio.run(server.run())
```

### 10.3 Claude Desktop設定

**~/.claude/mcp_config.json** (macOS/Linux)
**%APPDATA%\Claude\mcp_config.json** (Windows)

```json
{
  "mcpServers": {
    "pyscf-front": {
      "command": "python",
      "args": [
        "/path/to/pyscf-front/mcp_server/src/server.py"
      ],
      "env": {
        "PYTHONPATH": "/path/to/pyscf-front/backend/src"
      }
    }
  }
}
```

## 11. 統合開発環境セットアップ

### 11.1 VS Code設定

**.vscode/settings.json**
```json
{
  "python.defaultInterpreterPath": "${workspaceFolder}/venv/bin/python",
  "python.linting.enabled": true,
  "python.linting.pylintEnabled": false,
  "python.linting.flake8Enabled": true,
  "python.formatting.provider": "black",
  "python.testing.pytestEnabled": true,
  "python.testing.unittestEnabled": false,
  "dart.flutterSdkPath": "${env:HOME}/flutter",
  "[python]": {
    "editor.formatOnSave": true,
    "editor.codeActionsOnSave": {
      "source.organizeImports": "explicit"
    }
  },
  "[dart]": {
    "editor.formatOnSave": true,
    "editor.rulers": [80],
    "editor.selectionHighlight": false,
    "editor.suggest.snippetsPreventQuickSuggestions": false,
    "editor.suggestSelection": "first",
    "editor.tabCompletion": "onlySnippets",
    "editor.wordBasedSuggestions": "off"
  }
}
```

### 11.2 開発用Docker Compose

**docker-compose.dev.yml**
```yaml
version: '3.8'

services:
  # Pythonバックエンド
  backend:
    build:
      context: ./docker/backend
      dockerfile: Dockerfile.dev
    volumes:
      - ./backend:/app
      - ./protos:/protos
    environment:
      - PYTHONDONTWRITEBYTECODE=1
      - PYTHONUNBUFFERED=1
      - DATABASE_URL=mysql+pymysql://pyscf_user:pyscf_password@mysql:3306/pyscf_dev
      - GRPC_PORT=50051
      - MCP_SERVER_ENABLED=false
    ports:
      - "50051:50051"
      - "8000:8000"  # FastAPI
    depends_on:
      - mysql
    command: python -m uvicorn src.main:app --reload --host 0.0.0.0 --port 8000

  # MySQL
  mysql:
    image: mysql:8.0
    environment:
      MYSQL_ROOT_PASSWORD: root
      MYSQL_DATABASE: pyscf_dev
      MYSQL_USER: pyscf_user
      MYSQL_PASSWORD: pyscf_password
    ports:
      - "3306:3306"
    volumes:
      - mysql_data:/var/lib/mysql
      - ./database/schema:/docker-entrypoint-initdb.d

  # phpMyAdmin（開発用）
  phpmyadmin:
    image: phpmyadmin/phpmyadmin
    environment:
      PMA_HOST: mysql
      PMA_USER: root
      PMA_PASSWORD: root
    ports:
      - "8080:80"
    depends_on:
      - mysql

volumes:
  mysql_data:
```

### 11.3 Makefileの作成

**Makefile**
```makefile
.PHONY: help setup-backend setup-frontend setup test run-backend run-frontend docker-up docker-down clean

help:
	@echo "使用可能なコマンド:"
	@echo "  make setup          - 開発環境の完全セットアップ"
	@echo "  make setup-backend  - バックエンド環境のセットアップ"
	@echo "  make setup-frontend - フロントエンド環境のセットアップ"
	@echo "  make test          - 全テストの実行"
	@echo "  make run-backend   - バックエンドサーバーの起動"
	@echo "  make run-frontend  - フロントエンドの起動"
	@echo "  make docker-up     - Docker環境の起動"
	@echo "  make docker-down   - Docker環境の停止"
	@echo "  make clean         - キャッシュファイルの削除"

setup: setup-backend setup-frontend
	@echo "開発環境のセットアップが完了しました"

setup-backend:
	cd backend && python -m venv venv
	cd backend && ./venv/bin/pip install --upgrade pip
	cd backend && ./venv/bin/pip install -r requirements.txt
	cd backend && ./venv/bin/pip install -e .

setup-frontend:
	cd frontend && flutter pub get
	cd frontend && dart run build_runner build

test:
	cd backend && ./venv/bin/pytest
	cd frontend && flutter test

run-backend:
	cd backend && ./venv/bin/python -m uvicorn src.main:app --reload

run-frontend:
	cd frontend && flutter run -d macos  # または windows, linux

docker-up:
	docker-compose -f docker-compose.dev.yml up -d

docker-down:
	docker-compose -f docker-compose.dev.yml down

clean:
	find . -type f -name "*.pyc" -delete
	find . -type d -name "__pycache__" -delete
	find . -type d -name ".pytest_cache" -delete
	cd frontend && flutter clean
```

## 12. トラブルシューティング

### 12.1 Python関連

**問題**: Python 3.14でライブラリのインストールに失敗する
```bash
# 解決策1: ソースからビルド
pip install --no-binary :all: パッケージ名

# 解決策2: Python 3.13を使用
pyenv install 3.13.0
pyenv local 3.13.0
```

**問題**: PySCFのインポートエラー
```bash
# NumPy/SciPyの再インストール
pip uninstall numpy scipy
pip install --no-cache-dir numpy scipy
pip install --no-cache-dir pyscf
```

### 12.2 Flutter関連

**問題**: Flutter doctorでエラー
```bash
# Android toolchainエラーの場合
flutter doctor --android-licenses

# Xcodeエラーの場合（macOS）
sudo xcode-select --switch /Applications/Xcode.app/Contents/Developer
sudo xcodebuild -runFirstLaunch
```

### 12.3 gRPC関連

**問題**: protobufコンパイルエラー
```bash
# protoc-gen-dartが見つからない
export PATH="$PATH:$HOME/.pub-cache/bin"

# バージョン不整合
pip install grpcio-tools==1.73.0
```

### 12.4 GPU関連

**問題**: CUDAが認識されない
```bash
# 環境変数の設定
export CUDA_HOME=/usr/local/cuda
export PATH=$CUDA_HOME/bin:$PATH
export LD_LIBRARY_PATH=$CUDA_HOME/lib64:$LD_LIBRARY_PATH

# 確認
python -c "import torch; print(torch.cuda.is_available())"
```

### 12.5 データベース関連

**問題**: MySQL接続エラー
```sql
-- 認証プラグインの問題
ALTER USER 'pyscf_user'@'localhost' IDENTIFIED WITH mysql_native_password BY 'pyscf_password';
FLUSH PRIVILEGES;
```

---

## 付録A: 開発環境チェックリスト

- [ ] Python環境（3.13.0 または 3.14.0b2）
- [ ] Flutter SDK 3.32.0
- [ ] MySQL 8.0
- [ ] Docker Desktop
- [ ] VS Code + 拡張機能
- [ ] Git設定
- [ ] gRPCツール
- [ ] プロジェクト構造の作成
- [ ] 依存関係のインストール
- [ ] データベースの初期化
- [ ] サンプルコードの動作確認
- [ ] GPU環境（オプション）
- [ ] MCPサーバー（オプション）

## 付録B: 便利なコマンド集

```bash
# Python仮想環境の再作成
deactivate
rm -rf venv
python -m venv venv
source venv/bin/activate

# Flutterの完全リセット
flutter clean
flutter pub cache clean
flutter pub get

# Dockerの完全クリーンアップ
docker system prune -a --volumes

# gRPCコードの再生成
python -m grpc_tools.protoc -I./protos --python_out=./backend/src/api --grpc_python_out=./backend/src/api ./protos/*.proto

# データベースのバックアップ
mysqldump -u pyscf_user -p pyscf_dev > backup_$(date +%Y%m%d).sql
```

---

**文書情報**
- バージョン: 1.0
- 作成日: 2025-06-13
- 作成者: PySCF_Front開発チーム