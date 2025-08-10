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
- **パッケージ管理:** npm (Node.js), uv (Python)

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

3.  Pythonの依存関係をインストールします。（`uv` が必要です）
    ```bash
    cd src/python
    uv sync
    cd ../..
    ```

4.  開発モードでアプリケーションを起動します。
    ```bash
    npm run dev
    ```
    これにより、フロントエンドとバックエンドがホットリロード付きで起動します。

## 📦 アプリケーションのパッケージ化

プラットフォームに応じた配布用のアプリケーションをビルドするには、以下のコマンドを実行します。

```bash
npm run package
これにより release ディレクトリに実行可能ファイルが生成されます。