# PySCF Native App MCP Server

PySCF Native Appの機能をClaude Desktopから利用できるようにするMCP（Model Context Protocol）サーバーです。

## 機能

### 分子検索・変換
- **searchPubChem**: PubChemデータベースで化合物を検索し、3D構造（XYZ形式）を取得
- **convertSmiles**: SMILES文字列を3D分子構造に変換
- **validateXYZ**: XYZ形式データの妥当性検証

### 量子化学計算
- **startCalculation**: DFT、HF、MP2、CCSD、TDDFTなどの量子化学計算を開始
- **listCalculations**: 計算履歴の一覧表示
- **getCalculationDetails**: 特定計算の詳細結果取得
- **getOrbitals**: 分子軌道情報とエネルギーレベル取得
- **getOrbitalCube**: 軌道可視化用CUBEファイル生成
- **getIRSpectrum**: 理論振動スペクトラム生成

### システム管理
- **getSupportedParameters**: サポートされる計算手法・基底関数・汎関数一覧
- **getSettings**: アプリケーション設定取得
- **updateSettings**: 並列計算数・リソース制限設定
- **getResourceStatus**: CPU・メモリ使用状況監視
- **testConnection**: PySCF Native Appサーバーとの接続テスト

## セットアップ手順

### 1. 依存関係のインストール

```bash
cd mcp-server
npm install
```

### 2. TypeScript型定義の生成

```bash
npm run type-gen
```

### 3. ビルド

```bash
npm run build
```

### 4. Claude Desktop設定

最も簡単な設定方法（推奨）：

```bash
cd mcp-server
npm run setup-config
```

このスクリプトが以下を自動実行します：
- プロジェクトパスを自動検出
- Claude Desktop設定ファイルを生成
- コピー＆ペーストして使用する設定内容を表示

生成された設定をClaude Desktopの設定ファイル（`~/Library/Application Support/Claude/claude_desktop_config.json`）に追加し、Claude Desktopを再起動してください。

<details>
<summary>手動設定方法（上級者向け）</summary>

Claude Desktopの設定ファイルに以下の形式で追加:

```json
{
  "mcpServers": {
    "pyscf-native": {
      "command": "node",
      "args": ["{{PROJECT_PATH}}/mcp-server/dist/index.js"],
      "env": {}
    }
  }
}
```

**注意**: `{{PROJECT_PATH}}`部分を実際のプロジェクトパスに置換してください。

</details>

## 使用方法

### 1. PySCF Native Appの起動

まず、PySCF Native Appを起動します:

```bash
cd /Users/goodapple/workspace/PySCF_native_app
npm run dev
```

アプリが`http://127.0.0.1:5000`で起動することを確認してください。

### 2. Claude Desktop での利用

Claude Desktopを再起動すると、PySCF Native App関連のツールが利用可能になります。

#### 使用例

**化合物検索:**
```
水の分子構造をPubChemから検索してください
```

**SMILES変換:**
```
ベンゼンのSMILES「c1ccccc1」を3D構造に変換してください
```

**量子化学計算の開始:**
```
メタン分子のDFT計算を開始してください。基底関数は6-31G(d)、汎関数はB3LYPを使用してください
```

**計算結果の確認:**
```
実行中の計算の一覧を表示してください
```

**分子軌道の可視化:**
```
計算ID「calc_xxx」のHOMO軌道のCUBEファイルを生成してください
```

## 開発

### ファイル構造

```
mcp-server/
├── src/
│   ├── index.ts              # メインサーバー
│   ├── client.ts             # PySCF API クライアント
│   ├── types.ts              # OpenAPI 型定義（自動生成）
│   └── tools/                # MCPツール実装
│       ├── pubchem.ts        # 分子検索・変換ツール
│       ├── quantum.ts        # 量子化学計算ツール
│       ├── system.ts         # システム管理ツール
│       └── index.ts          # ツール統合
├── package.json
├── tsconfig.json
├── claude_desktop_config.json
└── README.md
```

### 開発用コマンド

```bash
# 型定義を再生成
npm run type-gen

# 開発モードで実行
npm run dev

# ビルド
npm run build

# 本番モードで実行
npm start
```

## トラブルシューティング

### 接続エラー

**症状**: "PySCF Native Appサーバーに接続できません"

**解決方法**:
1. PySCF Native Appが起動していることを確認
2. `http://127.0.0.1:5000`にアクセスして正常に動作することを確認
3. ファイアウォール設定を確認

### Claude Desktop でツールが表示されない

**解決方法**:
1. `claude_desktop_config.json` のパスが正しいことを確認
2. Claude Desktop を完全に再起動
3. ログを確認（通常は Claude Desktop のデベロッパーツールで確認可能）

### ビルドエラー

**解決方法**:
1. Node.js バージョンが18以上であることを確認
2. 依存関係を再インストール: `rm -rf node_modules package-lock.json && npm install`
3. 型定義を再生成: `npm run type-gen`

## ライセンス

Apache-2.0

## 貢献

バグ報告や機能要望は、メインプロジェクトのIssueでお知らせください。