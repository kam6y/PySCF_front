# my-review 最小仕様化 - 設計

## 方針
- 複雑な引数解析を廃止し、必要最小限の `--timeout` のみ解釈する。
- stdin の有無で `codex review` と `codex exec` を切り替える。
- タイムアウト処理は既存の 3段フォールバック（timeout / gtimeout / python3）を保持する。

## フロー
1. `MY_REVIEW_TIMEOUT` か `--timeout` から `timeout_seconds` を決定。
2. `stdin` が TTY かどうかでメモ有無を判定。
3. `stdin` が空:
   - `codex review --uncommitted` を実行。
4. `stdin` がある:
   - 一時ファイルにレビュー依頼プロンプトを生成。
   - `codex exec -` で stdin として渡す。
5. `timeout` / `gtimeout` / `python3` で期限付き実行。

## 例（擬似コード）
```bash
set -euo pipefail

timeout_seconds=${MY_REVIEW_TIMEOUT:-900}
parse --timeout
validate timeout_seconds

stdin_data=""
if [ ! -t 0 ]; then stdin_data=$(cat); fi

if [ -n "$stdin_data" ]; then
  prompt_file=$(mktemp)
  write prompt + git status/diff
  run_cmd=(codex exec -)
  run_with_timeout < prompt_file
else
  run_cmd=(codex review --uncommitted)
  run_with_timeout
fi
```

## 変更点
- `--commit` / `--title` / `--base` 等の分岐削除。
- 引数透過と複雑なプロンプト合成の廃止。
- ロジックを 1 ファイル内で直線的に読み下せる構造に整理。
