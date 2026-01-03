#!/usr/bin/env bash
set -euo pipefail

timeout_seconds="${MY_REVIEW_TIMEOUT:-900}"

while [ $# -gt 0 ]; do
  case "$1" in
    --timeout)
      shift
      if [ $# -eq 0 ]; then
        echo "missing value for --timeout" >&2
        exit 2
      fi
      timeout_seconds="$1"
      ;;
    --timeout=*)
      timeout_seconds="${1#*=}"
      ;;
    --)
      shift
      if [ $# -gt 0 ]; then
        echo "unexpected arguments: $*" >&2
        exit 2
      fi
      break
      ;;
    *)
      echo "unexpected argument: $1" >&2
      exit 2
      ;;
  esac
  shift
done

if ! [[ "$timeout_seconds" =~ ^[0-9]+$ ]] || [ "$timeout_seconds" -le 0 ]; then
  echo "invalid timeout: $timeout_seconds" >&2
  exit 2
fi

if ! command -v codex >/dev/null 2>&1; then
  echo "codex command not found in PATH" >&2
  exit 127
fi

stdin_data=""
if [ ! -t 0 ]; then
  stdin_data="$(cat)"
fi

prompt_file=""
run_cmd=()

if [ -n "$stdin_data" ]; then
  prompt_file="$(mktemp "${TMPDIR:-/tmp}/my-review.XXXXXX")"
  trap 'rm -f "$prompt_file"' EXIT

  {
    echo "未コミット差分のレビューをお願いします。"
    echo ""
    echo "ユーザーメモ:"
    printf '%s\n' "$stdin_data"
    echo ""
    echo "git status --porcelain=v1:"
    git status --porcelain=v1 || true
    echo ""
    echo "git diff (unstaged):"
    git diff || true
    echo ""
    echo "git diff --staged:"
    git diff --staged || true
  } > "$prompt_file"

  run_cmd=(codex exec -)
else
  run_cmd=(codex review --uncommitted)
fi

run_with_timeout() {
  local timeout_bin="$1"
  local status

  if [ -n "$prompt_file" ]; then
    if "$timeout_bin" "$timeout_seconds" "${run_cmd[@]}" < "$prompt_file"; then
      status=0
    else
      status=$?
    fi
  else
    if "$timeout_bin" "$timeout_seconds" "${run_cmd[@]}"; then
      status=0
    else
      status=$?
    fi
  fi

  if [ "$status" -eq 124 ]; then
    echo "my-review timed out after ${timeout_seconds} seconds. Retry with MY_REVIEW_TIMEOUT=1800 npm run my-review or npm run my-review -- --timeout 1800" >&2
  fi

  return "$status"
}

if command -v timeout >/dev/null 2>&1; then
  run_with_timeout timeout
elif command -v gtimeout >/dev/null 2>&1; then
  run_with_timeout gtimeout
else
  python3 - "$timeout_seconds" "${prompt_file:-}" "${run_cmd[@]}" <<'PY'
import subprocess
import sys

timeout_seconds = int(sys.argv[1])
prompt_path = sys.argv[2]
cmd = sys.argv[3:]
stdin = sys.stdin
if prompt_path:
  try:
    stdin = open(prompt_path, "r")
  except OSError:
    pass
try:
  proc = subprocess.Popen(cmd, stdin=stdin, stdout=sys.stdout, stderr=sys.stderr)
  proc.wait(timeout=timeout_seconds)
  sys.exit(proc.returncode)
except subprocess.TimeoutExpired:
  proc.kill()
  proc.wait()
  sys.stderr.write(
      f"my-review timed out after {timeout_seconds} seconds. "
      "Retry with MY_REVIEW_TIMEOUT=1800 npm run my-review or "
      "npm run my-review -- --timeout 1800\n"
  )
  sys.exit(124)
PY
fi
