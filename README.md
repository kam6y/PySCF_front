# PySCF_front

これは、**Electron**、**React (TypeScript)**、**Python (Flask)** を使用して構築された、量子化学計算アプリケーションです。

PySCFとRDKitをバックエンドに利用し、分子構造の可視化、PubChemからの分子データ取得、そしてDFT（密度汎関数理論）計算などを実行できるデスクトゥトップアプリケーションを目指して開発しています。



## 🌟 主な機能

- **3D分子可視化:** 3Dmol.jsを利用して、分子構造をインタラクティブに表示・操作できます。
- **分子構造の取得:**
  - PubChemのデータベースから、化合物の名称やCIDで検索し、3D構造を取得します。
  - SMILES形式の文字列を3D構造に変換します。
- **量子化学計算:**
  - PySCFを利用して、DFT計算などの量子化学計算を実行します。
  - 計算結果（HOMO/LUMO軌道、SCFエネルギーなど）を表示します。
- **計算履歴の管理:** 過去の計算結果を一覧表示し、名前の変更や削除が可能です。
- **自動環境構築:** ワンコマンドで開発環境をセットアップできます（`npm run setup-env`）。
- **環境検証機能:** Python依存関係と環境の健全性を自動チェックします（`npm run verify-env`）。
- **統一実行環境:** 開発・本番環境で同一のGunicornベースサーバーを使用し、環境差異問題を解決。
- **設定管理システム:** JSON設定ファイルでサーバー挙動を一元管理（`npm run debug:config`で確認可能）。
- **包括的テスト機能:** ビルド、パッケージング、サーバー起動の各段階を個別テスト可能。

## 🛠️ 技術スタック

- **フロントエンド:** React, TypeScript
- **バックエンド:** Python, Flask, Gunicorn, PySCF, RDKit
- **デスクトップフレームワーク:** Electron
- **ビルドツール:** Webpack, Electron Builder, PyInstaller
- **パッケージ管理:** npm (Node.js), conda (Python)
- **実行環境:** 統一されたGunicornベースサーバー（開発・本番両環境）
- **設定管理:** JSON設定ファイルによる一元管理

## 📸 アプリケーション画面

![PySCF Front アプリケーション画面](PySCF_front_view.png)

## 🚀 開発の始め方

### 🎯 クイックスタート（推奨）

最も簡単な方法で開発環境をセットアップできます：

1.  リポジトリをクローンします。
    ```bash
    git clone [https://github.com/kam6y/Pyscf_front.git](https://github.com/kam6y/Pyscf_front.git)
    cd Pyscf_front
    ```

2.  Node.jsの依存関係をインストールします。
    ```bash
    npm install
    ```

3.  **自動環境構築スクリプトを実行します。**
    ```bash
    npm run setup-env
    ```
    
    このコマンドは以下を自動的に実行します：
    - conda環境の存在確認
    - 必要に応じてconda環境の作成
    - すべての依存関係のインストール
    - 環境の検証

4.  **環境とサーバー設定を確認します（推奨）。**
    ```bash
    # 環境の健全性をチェック
    npm run verify-env
    
    # サーバー設定を確認
    npm run debug:config
    ```

5.  開発モードでアプリケーションを起動します。
    ```bash
    conda activate pyscf-env
    npm run dev
    ```

### 🔧 手動セットアップ

conda環境を手動で管理したい場合：

1.  リポジトリをクローンします。
    ```bash
    git clone [https://github.com/kam6y/Pyscf_front.git](https://github.com/kam6y/Pyscf_front.git)
    cd Pyscf_front
    ```

2.  Node.jsの依存関係をインストールします。
    ```bash
    npm install
    ```

3.  Python環境を設定します。（**conda環境が必須です**）
    ```bash
    # Miniforgeのインストール (Apple Silicon Macの場合)
    curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-arm64.sh"
    bash Miniforge3-MacOSX-arm64.sh -b -p $HOME/miniforge3
    
    # conda環境をenvironment.ymlから作成（全ての依存関係を含む）
    source $HOME/miniforge3/etc/profile.d/conda.sh
    conda env create -f .github/environment.yml
    
    # 環境をアクティブ化
    conda activate pyscf-env
    ```
    
    > **重要**: このプロジェクトはconda環境での開発が必須です。conda環境が正しく設定されていない場合、アプリケーションはエラーダイアログを表示します。

### conda環境のカスタムパス設定

アプリケーションは以下の順序でconda環境を自動検出します：

1. **環境変数** `CONDA_ENV_PATH` （最優先）
2. **conda info --base** を使用したベースパスからの検出
3. **conda info --envs** を使用した環境一覧からの検出
4. **fallback候補パス**: `~/miniforge3`、`~/miniconda3`、`~/anaconda3`、`~/mambaforge`

**カスタムパスを使用する場合:**
```bash
# conda環境が標準的でない場所にある場合
export CONDA_ENV_PATH="/path/to/your/pyscf-env"

# 例: カスタムインストール場所
export CONDA_ENV_PATH="/opt/miniconda3/envs/pyscf-env"
```

**トラブルシューティング:**
- conda環境が見つからない場合は、まず `conda info --envs` で `pyscf-env` が表示されることを確認してください
- 環境変数を設定した場合は、ターミナルを再起動してから `npm run dev` を実行してください
- conda がインストールされていない場合は、上記のMiniforgeインストール手順に従ってください

4.  **環境と設定の検証（推奨）**
    ```bash
    conda activate pyscf-env
    
    # 環境の健全性をチェック
    npm run verify-env
    
    # サーバー設定ファイルを確認
    npm run debug:config
    
    # Gunicornサーバーをローカルテスト
    npm run test:gunicorn-local
    ```
    
    環境に問題がある場合、詳細な診断情報とトラブルシューティング手順が表示されます。

5.  開発モードでアプリケーションを起動します。
    ```bash
    # conda環境をアクティブ化
    conda activate pyscf-env
    
    # アプリケーションの起動（統一Gunicornベースサーバーを使用）
    npm run dev
    ```
    これにより、フロントエンドとバックエンドが統一されたGunicornサーバーでホットリロード付きで起動します。

## 📦 アプリケーションのパッケージ化

### 🎁 配布用ビルド

プラットフォームに応じた配布用のアプリケーションをビルドするには、以下のコマンドを実行します。

```bash
# 完全ビルドテスト（推奨）
npm run test:build

# パッケージ作成
npm run package

# 完全なパッケージングテスト
npm run test:run-packaged
```

このコマンドは以下を自動的に実行します：
- フロントエンドのプロダクションビルド
- **統一Gunicorn環境の完全パッケージ化** (conda-packを使用)
- Python実行ファイルの生成 (PyInstallerを使用、Gunicorn対応)
- 統一サーバー設定ファイルの同梱
- Electronアプリケーションの配布パッケージ作成

### ✨ 配布の特徴

**エンドユーザーにとっての利点：**
- 🚫 **condaのインストール不要** - 完全な環境が同梱されています
- 🚫 **Python環境構築不要** - すべての依存関係が含まれています  
- ⚡ **即座に実行可能** - インストール後すぐに使用できます
- 🔒 **環境の隔離** - システムの Python 環境に影響しません
- 🎯 **統一実行環境** - 開発環境と同じGunicornベースサーバーで一貫性を保証
- ⚙️ **設定の一元管理** - サーバー挙動が設定ファイルで制御され、トラブルシューティングが容易

### 📂 生成されるファイル

```
dist/
├── PySCF_front-darwin-arm64.dmg     # macOS用
├── PySCF_front-win32-x64.exe        # Windows用  
└── PySCF_front-linux-x86_64.AppImage # Linux用
```

**内部構造：**
- `conda_env/` - 完全なPython環境（PySCF、RDKit、Gunicorn含む）
- `python_dist/` - PyInstaller実行ファイル（Gunicorn対応、フォールバック用）
- `config/` - 統一サーバー設定ファイル
- Electron アプリケーション

## 🛠️ トラブルシューティング

### 🔍 問題の診断

まず環境検証コマンドで問題を特定してください：

```bash
conda activate pyscf-env

# 環境の包括的チェック
npm run verify-env

# サーバー設定の確認
npm run debug:config

# 各コンポーネントの個別テスト
npm run test:python-build    # Python依存関係テスト
npm run test:gunicorn-local  # Gunicornサーバーテスト
npm run test:build          # 完全ビルドテスト
```

これらのコマンドは以下をチェックし、詳細な診断情報を提供します：
- Python バージョンと依存関係（Gunicorn含む）
- PySCF と RDKit の動作確認
- Flask と WebSocket 機能
- Gunicorn サーバーの動作確認
- サーバー設定ファイルの妥当性
- conda 環境の状態
- プロジェクト構造の整合性

### ⚡ よくある問題と解決方法

#### 1. 環境構築の問題

**問題:** `conda environment 'pyscf-env' not found`

**解決方法:**
```bash
# 自動環境構築を試行
npm run setup-env

# または手動で環境を作成
conda env create -f .github/environment.yml
conda activate pyscf-env
```

**問題:** `conda command not found`

**解決方法:**
```bash
# Miniforge をインストール
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-arm64.sh"
bash Miniforge3-MacOSX-arm64.sh -b -p $HOME/miniforge3
source $HOME/miniforge3/etc/profile.d/conda.sh
```

#### 2. 開発サーバーの問題

**問題:** `Python backend failed to start` または `Gunicorn startup failed`

**解決方法:**
1. 包括的環境の検証:
   ```bash
   conda activate pyscf-env
   npm run verify-env
   npm run debug:config
   ```
2. Gunicornサーバーの個別テスト:
   ```bash
   npm run test:gunicorn-local
   ```
3. 手動でのPythonサーバーテスト:
   ```bash
   conda activate pyscf-env
   cd src/python
   python app.py  # 直接実行テスト
   ```

**問題:** `Port already in use` や接続エラー

**解決方法:**
- アプリケーションは設定ファイルに基づいて自動的に空いているポートを検出します
- 設定ファイルの確認: `npm run debug:config`
- ファイアウォールの設定を確認してください
- アプリケーションを完全に終了してから再起動してください

**問題:** `"Works in dev but not in production"` または環境差異エラー

**解決方法:**
- **この問題は統一実行環境により解決されています**
- 開発・本番両環境で同じGunicornベースサーバーを使用
- 設定ファイルの妥当性確認: `npm run debug:config`
- 包括的テスト: `npm run test:build`

#### 3. パッケージ依存関係の問題

**問題:** `ModuleNotFoundError` や依存関係エラー

**解決方法:**
```bash
# conda環境を再作成
conda env remove -n pyscf-env
npm run setup-env

# または手動で依存関係を更新
conda activate pyscf-env
conda env update -f .github/environment.yml
```

#### 4. ビルドとパッケージ化の問題

**問題:** `PyInstaller build failed` または `Gunicorn dependency missing`

**解決方法:**
1. conda 環境をアクティベート: `conda activate pyscf-env`
2. 依存関係の確認:
   ```bash
   npm run test:python-build  # Python依存関係テスト
   npm run verify-env         # 環境の包括チェック
   ```
3. ビルドディレクトリをクリーン: `rimraf python_dist build/pyinstaller`
4. 再ビルド: `npm run build:python`

**問題:** パッケージ化後のアプリケーションが起動しない

**解決方法:**
- **統一実行環境テスト**:
  ```bash
  npm run test:build          # 完全ビルドテスト
  npm run test:run-packaged   # パッケージングテスト
  ```
- **診断:** パッケージには両方の実行方法が含まれています
  - アプリケーションはまず同梱conda環境（Gunicorn含む）を試行
  - 失敗時にPyInstaller実行ファイル（Gunicorn対応）にフォールバック
  - 統一サーバー設定により一貫した動作を保証

### 🆘 サポートが必要な場合

1. **環境情報の収集:**
   ```bash
   # 包括的環境検証の実行（詳細出力）
   npm run verify-env
   
   # サーバー設定の確認
   npm run debug:config
   
   # 各コンポーネントのテスト
   npm run test:python-build
   npm run test:gunicorn-local
   npm run test:build
   
   # システム情報
   conda info --envs
   python --version
   node --version
   npm --version
   ```

2. **ログの確認:**
   - アプリケーション起動時のコンソール出力
   - `build/pyinstaller/` ディレクトリのログファイル
   - Gunicornサーバーのログ出力
   - サーバー設定ファイルの内容

3. **問題報告:** 
   Issue報告時は上記の情報を含めてください: [GitHub Issues](https://github.com/kam6y/Pyscf_front/issues)