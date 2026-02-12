This file defines operating instructions for coding agents in this repository.

## 目的
- このプロジェクトは `PySCF_front`（Electron + React + Flask）です。
- 主機能は分子可視化、量子化学計算、計算結果管理、AIチャット支援です。
- 本ドキュメントの目的は、実装現状に一致した最小限かつ実用的な作業規約を提供することです。
- 仕様の一次情報は実装コードと `src/api-spec/openapi.yaml` です。古い文書よりコードを優先します。

## 適用範囲と優先順位
- `AGENTS.md` は実装・仕様に関する共通規約の正本です。
- `CLAUDE.md` は Claude Code の司令塔運用（委譲・レビュー手順）の正本です。
- 実装判断で矛盾した場合は、`実装コード` → `src/api-spec/openapi.yaml` → `AGENTS.md` の順で優先します。
- Claude専用の運用手順は `CLAUDE.md` を参照します。

## 会話ルール
- 返信・進捗報告・最終報告は常に日本語で行います。
- 変更前に「何を直すか」を短く共有し、変更後に結果と検証内容を明示します。
- 推測が含まれる場合は仮定として明示し、確認した事実と分けて記述します。
- 実行できなかった検証（時間不足、環境不足など）は必ず理由付きで報告します。

## 開発原則
- 開発段階のため、後方互換性より設計の単純性を優先します。
- breaking change は許容します。旧実装への過剰な互換レイヤーは追加しません。
- 複雑なフォールバックや使われていないコードは積極的に削除します。
- API変更時は、実装・型・ドキュメントを同一ターンで揃えます。
- ユーザーが明示的に求めない限り、破壊的なGit操作（`reset --hard` など）は行いません。

## ソースオブトゥルース
- API契約: `src/api-spec/openapi.yaml`
- Python生成モデル: `src/python/generated_models.py`（`npm run codegen` で再生成）
- TypeScript生成型: `src/web/types/generated-api.ts`（`npm run codegen` で再生成）
- サーバー設定: `config/server-config.json`
- 実行スクリプト: `package.json` の `scripts`
- AIチャット実装の現仕様: `src/python/api/agent.py`（Geminiベース単一チャット + SSE + chat history）

整合性ルール:
1. 実装と文書が矛盾したら、まず実装を確認し、文書を更新して一致させる。
2. OpenAPIに関わる変更は `npm run codegen` を前提に扱う。
3. 生成ファイルを手編集しない（必要なら生成元を修正する）。

## 変更フロー
1. 対象機能の実装箇所・型・API定義を確認する。
2. 変更設計を最小単位に分割し、影響範囲（UI/Backend/API/設定）を明確化する。
3. 実装を変更する。
4. 必要な生成・検証コマンドを実行する。
5. 変更差分を自己レビューし、古い説明・未使用コード・不要分岐を除去する。
6. 最終報告では「変更ファイル」「挙動への影響」「実行した検証」「未実施項目」を示す。

OpenAPI関連の標準順序:
1. `src/api-spec/openapi.yaml` を更新
2. `npm run codegen`
3. Python/Frontend実装を更新
4. テスト更新

## 最小コマンド
日常開発で使う最小セット（`package.json` と一致するもののみ記載）:

```bash
# セットアップ・検証
npm run setup-env
npm run verify-env
npm run verify-build-env

# 開発
npm run dev
npm run codegen
npm run debug:config

# ビルド・配布
npm run build
npm run package
npm run package:linux
npm run package:linux:docker
```

補助コマンド:
```bash
# 主要生成物を含むLinux向けビルド
npm run build:linux

# Pythonスタンドアロン確認
npm run test:python-standalone
```

## テスト実行規約
- Pythonテストは conda 環境の Python を優先して実行します。
- `python`/`python3` のシステム実行は依存差異を生みやすいため、原則避けます。

推奨手順:
```bash
conda activate pyscf-env
cd src/python
python -m pytest tests/ -v
```

直接パス実行の例:
```bash
~/miniforge3/envs/pyscf-env/bin/python -m pytest tests/ -v
~/miniforge3/envs/pyscf-env/bin/python -m pytest tests/integration/test_api_endpoints/test_quantum_api.py -v
```

実務ルール:
- 変更に近いテストを優先し、必要に応じて対象を絞る（`-k`, 単一ファイル, 単一テスト）。
- テスト未実行で提出する場合は、未実行理由と想定リスクを報告する。

## 外部ライブラリ調査規約
- 外部ライブラリを使う実装前に、最新APIを確認してからコードを書く。
- 特に更新頻度が高い領域（LLM周辺、フロントエンド、Flask拡張）は毎回確認する。

必須確認項目:
1. APIシグネチャ（引数名・戻り値・型）
2. 非推奨/削除済みAPI
3. 現行の推奨実装パターン
4. 直近バージョンの破壊的変更

調査方針:
- 一次情報（公式ドキュメント、公式リポジトリ、公式リリースノート、context7）を優先する。
- 参考記事は補助扱いにし、公式情報で再確認する。
- 実装時は参照元を簡潔に記録し、なぜその書き方を採用したか説明可能にする。

## ドメイン注意点
- 量子化学計算は `DFT`, `HF`, `MP2`, `CCSD`, `CCSD_T`, `TDDFT`, `CASCI`, `CASSCF` を扱う。
- リアルタイム更新は Socket.IO のイベント駆動（`join_calculation`, `leave_calculation`, `join_global_updates`, `calculation_update`）で行う。
- AI機能の現仕様は「Gemini APIベースの単一チャット」です。
- `/api/agent/chat` は SSE で `agent_status` / `chunk` / `done` / `error` を返す。
- チャット履歴は `session_id` がある場合にDBへ保存する。
- GPU4PySCF は Linux のみ対象。利用不可・失敗時は CPU にフォールバックする。
- PySCFの `spin` は不対電子数（2S）であり、一般的な多重度（2S+1）と表記が異なる。

## 更新時チェックリスト
- [ ] 返信を日本語で統一した。
- [ ] 変更内容が「単純化優先・breaking change許容」の方針に沿っている。
- [ ] API変更時に `src/api-spec/openapi.yaml` と生成物の整合を取った。
- [ ] 記載した `npm run ...` が `package.json` の `scripts` に存在する。
- [ ] GPU説明が Linux 限定かつ CPUフォールバック前提になっている。
- [ ] 実行した検証コマンドと未実施項目を最終報告に記載した。
