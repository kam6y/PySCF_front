## 1. プロジェクト概要

### 1.1 背景と目的
PySCF_Frontは、PySCF（Python-based Simulations of Chemistry Framework）に対する使いやすいGUIフロントエンドを提供するクロスプラットフォームアプリケーションです。

量子化学計算は研究・教育・産業分野で重要な役割を果たしているが、既存のソ# PySCF_Front 要件定義書
## 量子化学計算クロスプラットフォームアプリケーション

## 1. プロジェクト概要

### 1.1 背景と目的
PySCF_Frontは、PySCF（Python-based Simulations of Chemistry Framework）に対する使いやすいGUIフロントエンドを提供するクロスプラットフォームアプリケーションです。

量子化学計算は研究・教育・産業分野で重要な役割を果たしているが、既存のソフトウェアは以下の課題を抱えている：
- コマンドライン操作が中心で初学者にとって敷居が高い
- 複数の計算ジョブ管理が煩雑
- 可視化機能が別ツールに依存
- プラットフォーム間での互換性が低い

本プロジェクトは、これらの課題を解決し、直感的なGUIと強力な計算エンジンを組み合わせた統合環境を提供することを目的とする。

### 1.2 対象ユーザー
- 量子化学研究者
- 化学系の大学生・大学院生
- 創薬・材料開発に従事する企業研究者
- 計算化学の教育者

### 1.3 開発スコープ
- PySCFを基盤とした量子化学計算エンジンのGUIフロントエンド
- Electronによるクロスプラットフォーム対応デスクトップアプリケーション
- 3Dpymolによる分子可視化機能
- 拡張性を考慮したモジュラー設計
- 将来的なMCPサーバー実装によるClaude Desktop連携

## 2. 機能要件

### 2.1 分子入力機能
- **RDKit GUIによる分子描画**
  - 2D分子エディタによる直感的な構造入力
  - 自動的な3D構造生成と最適化
- **SMILES文字列入力**
  - テキストベースの分子記述
  - バッチ入力対応（複数SMILES一括処理）
- **PubChem連携**
  - 化合物名による検索
  - PubChem CIDによる直接指定
- **座標ファイル入力**
  - XYZ形式ファイルのインポート
  - MOL/SDF形式対応
  - PDB形式対応（タンパク質等の大分子用）

### 2.2 量子化学計算機能

#### 2.2.1 計算手法（優先順位順）
1. **DFT (Density Functional Theory)**
   - 汎関数：B3LYP, PBE, PBE0, M06-2X等
   - グリッド設定のカスタマイズ
2. **HF (Hartree-Fock)**
   - RHF, UHF, ROHF対応
3. **MP2 (Møller-Plesset 2次摂動)**
   - RI-MP2対応
4. **その他の手法**
   - CCSD(T)
   - CASSCF/CASPT2
   - TD-DFT（励起状態計算）

#### 2.2.2 計算タイプ
1. **1点エネルギー計算**（最優先実装）
2. **構造最適化**
   - 局所最小構造探索
   - 収束条件のカスタマイズ
3. **振動解析**
   - 振動数計算
   - IR/Ramanスペクトル予測
   - 熱力学的性質計算
4. **遷移状態探索**
5. **反応経路探索**
6. **分子動力学計算**

#### 2.2.3 基底関数セット
- Pople系：6-31G*, 6-311++G(d,p)等
- Dunning系：cc-pVDZ, aug-cc-pVTZ等
- Ahlrichs系：def2-SVP, def2-TZVP等
- カスタム基底関数の定義機能

### 2.3 計算ジョブ管理機能

#### 2.3.1 バッチ計算システム
- **ジョブキュー管理**
  - 計算待機リストの表示
  - 優先度設定
  - 一時停止/再開機能
- **並列計算制御**
  - CPU使用率の設定
  - メモリ使用量の制限
- **自動実行**
  - 前の計算終了を検知して次の計算を開始
  - エラー時の自動リトライ機能

#### 2.3.2 計算パラメータ管理
- **推奨設定テンプレート**
  - 分子サイズに応じた自動設定
  - 計算目的別プリセット
- **詳細設定**
  - SCF収束条件
  - 積分精度
  - 対称性の利用設定

### 2.4 可視化機能

#### 2.4.1 分子構造表示（3Dpymol使用）
- Ball-and-stick、Space-filling等の表示モード
- 原子ラベル、結合長・角度の表示
- アニメーション（振動モード、反応経路）

#### 2.4.2 計算結果の可視化
- 分子軌道（HOMO/LUMO等）
- 電子密度分布
- 静電ポテンシャルマップ
- スペクトル表示（IR、UV-Vis等）

### 2.5 プロジェクト管理機能
- **プロジェクト作成・保存**
  - 分子構造、計算設定、結果の一括管理
  - プロジェクトテンプレート機能
- **履歴管理**
  - 計算履歴の記録
  - 結果の比較機能
- **タグ・検索機能**
  - プロジェクトへのタグ付け
  - 分子名、計算手法等での検索

### 2.6 データ入出力機能

#### 2.6.1 内部データ管理
- **JSON形式での保存**
  - 分子構造情報
  - 計算パラメータ
  - 計算結果
  - PySCFログ全文

#### 2.6.2 エクスポート機能
- **ログファイル出力**（.qclog等の独自拡張子）
- **標準フォーマット出力**
  - Gaussian形式
  - ORCA形式
  - Molden形式
- **レポート生成**
  - PDF/HTML形式
  - 計算サマリー自動生成

## 3. 非機能要件

### 3.1 性能要件
- **最大原子数**：1000原子
- **計算速度**：PySCFネイティブ性能の90%以上
- **UI応答性**：計算中もUIが応答（非ブロッキング）

### 3.2 ユーザビリティ要件
- **直感的なUI**
  - ドラッグ&ドロップ対応
  - ツールチップによるヘルプ表示
- **多言語対応**
  - 日本語/英語切り替え
- **ダークモード対応**

### 3.3 信頼性要件
- **自動保存機能**
  - 5分ごとの自動保存
  - クラッシュリカバリ機能
- **エラーハンドリング**
  - 詳細なエラーメッセージ
  - エラーログの記録

### 3.4 互換性要件
- **対応OS**
  - Windows 10/11
  - macOS 10.15以降
  - Ubuntu 20.04 LTS以降
- **Python環境**
  - Python 3.13対応
  - 仮想環境の自動構築

## 4. システムアーキテクチャ

### 4.1 全体構成
```
┌─────────────────────────────────────────────┐
│          Electron Frontend                   │
│  ┌─────────────┐  ┌────────────────────┐   │
│  │   React UI  │  │  3Dpymol Viewer    │   │
│  └─────────────┘  └────────────────────┘   │
└───────────────────┬─────────────────────────┘
                    │ IPC (Context Bridge)
┌───────────────────┴─────────────────────────┐
│          Electron Main Process               │
│  ┌─────────────┐  ┌────────────────────┐   │
│  │ Job Manager │  │  Data Storage      │   │
│  └─────────────┘  └────────────────────┘   │
└───────────────────┬─────────────────────────┘
                    │ Child Process
┌───────────────────┴─────────────────────────┐
│          Python Backend (PySCF)              │
│  ┌─────────────┐  ┌────────────────────┐   │
│  │ Calculation │  │  Result Parser     │   │
│  │   Engine    │  │                    │   │
│  └─────────────┘  └────────────────────┘   │
└─────────────────────────────────────────────┘
```

### 4.3 PySCFとの統合

#### 4.3.1 PySCFの主要機能活用
- **基底関数セット**：PySCF内蔵の基底関数（STO-3G、6-31G*、cc-pVDZ、def2系など）を直接利用
- **交換相関汎関数**：Libxc（デフォルト）またはXCFunライブラリを通じて500以上の汎関数にアクセス
- **高度な機能**：
  - 密度フィッティング（DF/RI）による高速化
  - X2Cによるスカラー相対論効果の考慮
  - 自動微分機能（PySCFAD）の活用検討
  - GPU加速（GPU4PySCF）の将来的な統合

#### 4.3.2 計算フロー
1. **入力処理**：Electron側からの分子構造・計算条件をPySCFのMoleオブジェクトに変換
2. **計算実行**：PySCFの各種メソッド（scf.RHF、dft.RKS等）を呼び出し
3. **結果処理**：PySCFの出力を解析し、Electron側で表示可能な形式に変換
4. **エラー処理**：PySCFの例外を適切にキャッチし、ユーザーフレンドリーなメッセージに変換

## 5. フォルダ構成

```
pyscf-front/
├── src/
│   ├── main/                    # Electronメインプロセス
│   │   ├── index.ts            # エントリーポイント
│   │   ├── ipc/                # IPC通信ハンドラ
│   │   ├── services/           # ビジネスロジック
│   │   └── preload.ts          # contextBridge設定
│   │
│   ├── renderer/               # Electronレンダラープロセス
│   │   ├── App.tsx            # Reactアプリケーション
│   │   ├── components/         # UIコンポーネント
│   │   ├── hooks/             # カスタムフック
│   │   ├── stores/            # 状態管理
│   │   └── utils/             # ユーティリティ
│   │
│   └── python/                 # Pythonバックエンド
│       ├── main.py            # エントリーポイント
│       ├── calculations/       # 計算エンジン
│       │   ├── engine.py      # PySCF計算実行エンジン
│       │   ├── methods/       # 計算手法ラッパー
│       │   │   ├── dft.py     # DFT計算ラッパー
│       │   │   ├── hf.py      # HF計算ラッパー
│       │   │   ├── mp2.py     # MP2計算ラッパー
│       │   │   └── post_hf.py # その他のpost-HF手法
│       │   ├── functionals/   # DFT汎関数設定
│       │   │   ├── xc_functionals.py  # 交換相関汎関数定義
│       │   │   ├── custom_xc.py       # カスタム汎関数
│       │   │   └── functional_db.json # 汎関数データベース
│       │   └── basis/         # 基底関数設定
│       │       ├── basis_manager.py   # 基底関数管理
│       │       └── custom_basis.py    # カスタム基底関数
│       ├── parsers/           # 結果パーサー
│       │   ├── output_parser.py       # PySCF出力解析
│       │   ├── property_extractor.py  # 物性値抽出
│       │   └── format_converter.py    # フォーマット変換
│       ├── visualization/     # 3Dpymol連携
│       │   ├── mol_viewer.py          # 分子構造可視化
│       │   ├── orbital_viewer.py      # 軌道可視化
│       │   └── density_viewer.py      # 電子密度可視化
│       └── utils/             # ユーティリティ
│           ├── molecule_builder.py    # 分子構造構築
│           ├── job_manager.py         # ジョブ管理
│           └── resource_monitor.py    # リソース監視
│
├── assets/                     # 静的リソース
├── config/                     # 設定ファイル
│   ├── default_settings.json  # デフォルト計算設定
│   ├── xc_functionals.json    # 汎関数プリセット定義
│   └── calculation_presets.json # 計算プリセット
├── data/                       # データ保存用
│   ├── projects/              # プロジェクトデータ
│   ├── templates/             # テンプレート
│   └── cache/                 # キャッシュ
├── tests/                      # テストコード
│   ├── unit/                  # 単体テスト
│   ├── integration/           # 統合テスト
│   └── benchmarks/            # ベンチマーク
└── docs/                       # ドキュメント
    ├── user_guide/            # ユーザーガイド
    ├── api_reference/         # APIリファレンス
    └── development/           # 開発者向けドキュメント
```

### 5.1 主要ファイルの役割

#### Python側の重要ファイル

**functionals/xc_functionals.py**
```python
# DFT交換相関汎関数の管理・設定
class XCFunctionalManager:
    """
    PySCFのLibxc/XCFunとのインターフェース
    - 汎関数の選択と設定
    - カスタム汎関数の定義
    - ハイブリッド汎関数のパラメータ管理
    """
    def get_functional_info(self, xc_name):
        """汎関数の詳細情報を取得"""
        pass
    
    def set_custom_functional(self, definition):
        """カスタム汎関数の定義"""
        pass
```

**functionals/functional_db.json**
```json
{
  "functionals": {
    "B3LYP": {
      "type": "hybrid",
      "description": "Becke 3-parameter Lee-Yang-Parr",
      "definition": "HF*0.2 + .08*LDA + .72*B88, .81*LYP + .19*VWN",
      "category": "popular"
    },
    "PBE": {
      "type": "GGA",
      "description": "Perdew-Burke-Ernzerhof",
      "category": "popular"
    }
  },
  "custom_presets": {
    "organic_molecules": {
      "functional": "B3LYP",
      "basis": "6-31G*",
      "reason": "バランスの取れた精度と計算コスト"
    }
  }
}
```

**basis/basis_manager.py**
```python
# PySCFの基底関数管理
class BasisSetManager:
    """
    PySCFの基底関数セットとのインターフェース
    - 内蔵基底関数の利用
    - カスタム基底関数の定義
    - Basis Set Exchangeとの連携
    """
    def get_available_basis(self):
        """利用可能な基底関数リストを取得"""
        pass
```

### 6.1 プロジェクトデータ構造
```json
{
  "projectId": "uuid",
  "projectName": "Project Name",
  "createdAt": "2025-01-23T10:00:00Z",
  "molecules": [
    {
      "moleculeId": "uuid",
      "name": "Benzene",
      "inputType": "SMILES",
      "inputData": "c1ccccc1",
      "coordinates": [...],
      "calculations": [
        {
          "calculationId": "uuid",
          "method": "DFT",
          "functional": "B3LYP",
          "basisSet": "6-31G*",
          "status": "completed",
          "results": {...}
        }
      ]
    }
  ]
}
```

### 6.2 計算ジョブキュー構造
```json
{
  "queue": [
    {
      "jobId": "uuid",
      "moleculeId": "uuid",
      "calculationType": "optimization",
      "priority": 1,
      "status": "waiting",
      "createdAt": "2025-01-23T10:00:00Z",
      "parameters": {...}
    }
  ]
}
```

## 7. UI/UXデザイン指針

### 7.1 画面構成
- **メイン画面**：3ペイン構成
  - 左：プロジェクト/分子リスト
  - 中央：3D分子ビューア
  - 右：計算設定/結果表示
- **ジョブ管理画面**：計算キューの可視化
- **設定画面**：計算パラメータのカスタマイズ

### 7.2 操作フロー
1. プロジェクト作成/選択
2. 分子入力（複数方法から選択）
3. 計算手法・パラメータ設定
4. 計算実行（キューに追加）
5. 結果確認・可視化
6. エクスポート

## 8. 開発ロードマップ

### Phase 1: 基盤構築（2ヶ月）
- Electron + React環境構築
- Python連携基盤実装
- 基本的なUI実装

### Phase 2: コア機能実装（3ヶ月）
- PySCF統合
- 1点計算機能
- 分子入力機能（SMILES、XYZ）
- 基本的な可視化

### Phase 3: 拡張機能実装（3ヶ月）
- 構造最適化、振動解析
- バッチ計算システム
- プロジェクト管理機能
- 高度な可視化機能

### Phase 4: 最適化・安定化（2ヶ月）
- パフォーマンス最適化
- バグ修正
- ドキュメント整備

### Phase 5: 将来拡張（継続的）
- MCPサーバー実装
- Claude Desktop連携
- コミュニティプラグイン対応

## 9. 将来構想

### 9.1 MCPサーバー連携
- **Claude Desktopからの操作**
  - 自然言語による分子指定
  - 計算条件の自動設定
  - 結果の要約生成
- **実装イメージ**
  ```json
  {
    "mcpServers": {
      "pyscf-front": {
        "command": "node",
        "args": ["path/to/mcp-server.js"],
        "tools": [
          "calculate_molecule",
          "get_calculation_status",
          "analyze_results"
        ]
      }
    }
  }
  ```

### 9.2 機能拡張の可能性
- 機械学習による反応予測
- クラウド計算リソースとの連携
- 実験データとの統合解析
- VR/ARによる分子操作

## 10. セキュリティ・プライバシー考慮事項

### 10.1 データ保護
- ローカルストレージのみ使用（クラウド非依存）
- プロジェクトファイルの暗号化オプション

### 10.2 計算リソース管理
- CPU/メモリ使用量の制限機能
- 悪意のある入力データのバリデーション

## 11. テスト戦略

### 11.1 単体テスト
- 計算エンジンの正確性検証
- UIコンポーネントテスト

### 11.2 統合テスト
- Electron-Python間通信テスト
- エンドツーエンドシナリオテスト

### 11.3 性能テスト
- 大規模分子での計算性能
- メモリリーク検証

## 12. ドキュメント計画

### 12.1 ユーザー向け
- クイックスタートガイド
- 計算手法の解説
- トラブルシューティング

### 12.2 開発者向け
- APIリファレンス
- プラグイン開発ガイド
- コントリビューションガイドライン