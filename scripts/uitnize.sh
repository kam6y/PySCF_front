#!/bin/bash

# 引数があればそのディレクトリ、なければカレントディレクトリ
TARGET_DIR="${1:-.}"
OUTPUT_FILE="code.txt"

# tree用の除外パターン（表示のみに使用）
TREE_IGNORE="node_modules|dist|build|out|release|python_dist|__pycache__|.pytest_cache|.vscode|.idea|.claude|.serena|.DS_Store|*.log|conda_env|.git"

{
  printf "\n"
  printf "=================================================\n"
  printf "ディレクトリ構造\n"
  printf "=================================================\n\n"

  # treeコマンドが利用可能かチェック（表示のみに使用）
  if command -v tree &> /dev/null; then
    tree -I "$TREE_IGNORE" "$TARGET_DIR"
  else
    echo "※ treeコマンドが見つかりません。簡易表示を使用します。"
    echo ""
    # Git環境ならgit ls-filesを使用、そうでなければfindにフォールバック
    if git rev-parse --git-dir > /dev/null 2>&1; then
      {
        git ls-files "$TARGET_DIR"
        git ls-files --others --exclude-standard "$TARGET_DIR"
      } | sort -u | sed 's|/[^/]*$||' | sort -u | grep -v "^$" || echo "$TARGET_DIR"
    else
      echo "※ 警告: Gitリポジトリが検出されませんでした。findコマンドを使用します。" >&2
      find "$TARGET_DIR" -type d \
        ! -path "*/node_modules/*" \
        ! -path "*/dist/*" \
        ! -path "*/build/*" \
        ! -path "*/out/*" \
        ! -path "*/release/*" \
        ! -path "*/python_dist/*" \
        ! -path "*/__pycache__/*" \
        ! -path "*/.pytest_cache/*" \
        ! -path "*/.vscode/*" \
        ! -path "*/.idea/*" \
        ! -path "*/.claude/*" \
        ! -path "*/.serena/*" \
        ! -path "*/.git/*" \
        ! -path "*/conda_env/*" \
        | sort
    fi
  fi

  printf "\n\n"
  printf "=================================================\n"
  printf "ファイル内容\n"
  printf "=================================================\n\n"

  # ファイルを検索（.gitignoreに基づいて自動除外）
  if git rev-parse --git-dir > /dev/null 2>&1; then
    # Git環境: 追跡中のファイル + 未追跡だが.gitignoreで除外されないファイル
    files=$(
      {
        git ls-files "$TARGET_DIR"
        git ls-files --others --exclude-standard "$TARGET_DIR"
      } | sort -u | grep -v "^${OUTPUT_FILE}$"
    )
  else
    # Git非環境: findコマンドにフォールバック
    echo "※ 警告: Gitリポジトリが検出されませんでした。findコマンドを使用します。" >&2
    files=$(find "$TARGET_DIR" -type f \
      ! -path "*/node_modules/*" \
      ! -path "*/dist/*" \
      ! -path "*/build/*" \
      ! -path "*/out/*" \
      ! -path "*/release/*" \
      ! -path "*/python_dist/*" \
      ! -path "*/__pycache__/*" \
      ! -path "*/.pytest_cache/*" \
      ! -path "*/.vscode/*" \
      ! -path "*/.idea/*" \
      ! -path "*/.claude/*" \
      ! -path "*/.serena/*" \
      ! -path "*/.git/*" \
      ! -path "*/conda_env/*" \
      ! -name ".DS_Store" \
      ! -name "*.log" \
      ! -name "*.pyc" \
      ! -name "chat_history.db" \
      ! -name "code.txt" \
      | sort)
  fi

  for f in $files; do
    printf "%s\n" "----------------------------------------------------"
    printf "%s\n" "$f"
    printf "%s\n\n" "----------------------------------------------------"

    # バイナリファイルのチェック（拡張子ベース + fileコマンド）
    file_ext="${f##*.}"
    if [[ "$file_ext" =~ ^(txt|md|json|yaml|yml|js|jsx|ts|tsx|py|sh|css|html|xml|svg|csv|log|env|gitignore|prettierrc|prettierignore)$ ]] || file "$f" | grep -qi "text\|json\|xml\|script"; then
      nl -ba "$f"
    else
      echo "  [バイナリファイル - スキップ]"
    fi

    printf "\n\n"
  done
} > "$OUTPUT_FILE"

echo "✓ 出力が完了しました: $OUTPUT_FILE"
