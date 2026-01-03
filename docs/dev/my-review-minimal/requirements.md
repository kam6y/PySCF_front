# my-review 最小仕様化 - 要件

## 背景/目的
- `scripts/my-review.sh` が多機能化しているため、運用に必要な最小仕様に絞って簡素化する。

## 目標
- 必須用途（未コミット差分レビュー + stdinメモ）だけを確実に動作させる。
- 使い方は AGENTS.md に記載の例をそのまま満たす。

## 対象範囲
- `scripts/my-review.sh` のロジック簡素化
- `package.json` の `my-review` スクリプトは維持

## 非対象
- `--commit` / `--title` / `--base` などの詳細オプション
- 追加の codex オプション透過
- 高度な引数解析

## 機能要件
1. **デフォルト動作**
   - `stdin` が空の場合、`codex review --uncommitted` を実行する。
2. **stdin 指定時の動作**
   - `stdin` にメモが流れている場合、`codex exec -` を実行し、以下を含むレビュー依頼プロンプトを生成する。
     - 「未コミット差分のレビューをお願いします。」
     - 「ユーザーメモ:」+ stdin 内容
     - `git status --porcelain=v1`
     - `git diff` (unstaged)
     - `git diff --staged`
3. **タイムアウト**
   - 既定は `MY_REVIEW_TIMEOUT` 環境変数（未設定時は 900 秒）。
   - `--timeout` / `--timeout=SECONDS` を受け付ける。
4. **タイムアウト実装**
   - `timeout` → `gtimeout` → `python3` の順で利用する。
5. **エラーハンドリング**
   - `codex` コマンドが見つからない場合はエラー終了。
   - `--timeout` が不正（数値以外/0以下）の場合はエラー終了。

## 受け入れ基準
- `npm run my-review` で未コミット差分レビューが走る。
- `cat <<'EOF' | npm run my-review ... EOF` でメモ付きレビューが走る。
- `MY_REVIEW_TIMEOUT=1800 npm run my-review` が有効。
- `npm run my-review -- --timeout 1800` が有効。
- 既存の `scripts/my-review.sh` より短く読みやすい。
