set shell := ["bash", "-cu"]
# Codexサンドボックスでのrename制限回避のため、書き込み可能な固定パスを使用する
uv_cache_dir := "/home/grace/.codex/uv-cache"
uv_tmp_dir := "/home/grace/.codex/uv-tmp"
web_port := "3001"
api_port := "8000"

@default:
  just --list

web:
  #!/usr/bin/env bash
  set -euo pipefail
  port="${WEB_PORT:-{{web_port}}}"
  initial_port="$port"
  port_in_use() {
    if command -v ss >/dev/null 2>&1; then
      ss -ltn "sport = :$1" | awk 'NR>1 {found=1} END {exit !found}'
    elif command -v lsof >/dev/null 2>&1; then
      lsof -nP -iTCP:"$1" -sTCP:LISTEN >/dev/null 2>&1
    else
      return 1
    fi
  }
  find_free_port() {
    local candidate="$1"
    for _ in $(seq 0 20); do
      if ! port_in_use "$candidate"; then
        echo "$candidate"
        return 0
      fi
      candidate=$((candidate + 1))
    done
    echo "利用可能なポートが見つかりません (start=$1)" >&2
    return 1
  }
  port="$(find_free_port "$port")"
  if [[ "$port" != "$initial_port" ]]; then
    echo "WEB_PORT $initial_port は使用中のため $port を利用します"
  fi
  pnpm -C apps/web dev --port "$port"

api:
  #!/usr/bin/env bash
  set -euo pipefail
  port="${API_PORT:-{{api_port}}}"
  initial_port="$port"
  port_in_use() {
    if command -v ss >/dev/null 2>&1; then
      ss -ltn "sport = :$1" | awk 'NR>1 {found=1} END {exit !found}'
    elif command -v lsof >/dev/null 2>&1; then
      lsof -nP -iTCP:"$1" -sTCP:LISTEN >/dev/null 2>&1
    else
      return 1
    fi
  }
  find_free_port() {
    local candidate="$1"
    for _ in $(seq 0 20); do
      if ! port_in_use "$candidate"; then
        echo "$candidate"
        return 0
      fi
      candidate=$((candidate + 1))
    done
    echo "利用可能なポートが見つかりません (start=$1)" >&2
    return 1
  }
  port="$(find_free_port "$port")"
  if [[ "$port" != "$initial_port" ]]; then
    echo "API_PORT $initial_port は使用中のため $port を利用します"
  fi
  cd apps/api
  mkdir -p {{uv_cache_dir}} {{uv_tmp_dir}}
  TMPDIR={{uv_tmp_dir}} UV_CACHE_DIR={{uv_cache_dir}} uv sync
  TMPDIR={{uv_tmp_dir}} UV_CACHE_DIR={{uv_cache_dir}} uv run uvicorn main:app --reload --port "$port"

dev:
  #!/usr/bin/env bash
  set -euo pipefail
  trap 'kill 0' INT TERM EXIT
  port_in_use() {
    if command -v ss >/dev/null 2>&1; then
      ss -ltn "sport = :$1" | awk 'NR>1 {found=1} END {exit !found}'
    elif command -v lsof >/dev/null 2>&1; then
      lsof -nP -iTCP:"$1" -sTCP:LISTEN >/dev/null 2>&1
    else
      return 1
    fi
  }
  find_free_port() {
    local candidate="$1"
    for _ in $(seq 0 20); do
      if ! port_in_use "$candidate"; then
        echo "$candidate"
        return 0
      fi
      candidate=$((candidate + 1))
    done
    echo "利用可能なポートが見つかりません (start=$1)" >&2
    return 1
  }
  api_port="${API_PORT:-{{api_port}}}"
  web_port="${WEB_PORT:-{{web_port}}}"
  initial_api_port="$api_port"
  initial_web_port="$web_port"
  api_port="$(find_free_port "$api_port")"
  web_port="$(find_free_port "$web_port")"
  if [[ "$api_port" != "$initial_api_port" ]]; then
    echo "API_PORT $initial_api_port は使用中のため $api_port を利用します"
  fi
  if [[ "$web_port" != "$initial_web_port" ]]; then
    echo "WEB_PORT $initial_web_port は使用中のため $web_port を利用します"
  fi
  pushd apps/api >/dev/null
  mkdir -p {{uv_cache_dir}} {{uv_tmp_dir}}
  TMPDIR={{uv_tmp_dir}} UV_CACHE_DIR={{uv_cache_dir}} uv sync
  TMPDIR={{uv_tmp_dir}} UV_CACHE_DIR={{uv_cache_dir}} uv run uvicorn main:app --reload --port "$api_port" &
  popd >/dev/null
  pnpm -C apps/web dev --port "$web_port" &
  wait
