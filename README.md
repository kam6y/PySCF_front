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

## 🛠️ 技術スタック

- **フロントエンド:** React, TypeScript
- **バックエンド:** Python, Flask, PySCF, RDKit
- **デスクトップフレームワーク:** Electron
- **ビルドツール:** Webpack, Electron Builder, PyInstaller
- **パッケージ管理:** npm (Node.js), conda (Python)

## 🚀 開発の始め方

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
    
    # conda環境の作成とアクティブ化
    source $HOME/miniforge3/etc/profile.d/conda.sh
    conda create -y -n pyscf-env python=3.12
    conda activate pyscf-env
    
    # パッケージのインストール
    conda install -y -c conda-forge pyscf rdkit flask geometric requests flask-cors pydantic gevent threadpoolctl
    pip install flask-sock flask-pydantic datamodel-code-generator pyinstaller gevent-websocket certifi
    ```
    
    > **重要**: このプロジェクトはconda環境での開発が必須です。conda環境が正しく設定されていない場合、アプリケーションはエラーダイアログを表示します。

4.  開発モードでアプリケーションを起動します。
    ```bash
    # conda環境をアクティブ化
    conda activate pyscf-env
    
    # アプリケーションの起動
    npm run dev
    ```
    これにより、フロントエンドとバックエンドがホットリロード付きで起動します。

## 📦 アプリケーションのパッケージ化

プラットフォームに応じた配布用のアプリケーションをビルドするには、以下のコマンドを実行します。

```bash
npm run package
これにより PySCF_native_app/dist ディレクトリに実行可能ファイルが生成されます。
```