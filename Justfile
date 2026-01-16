set shell := ["bash", "-cu"]
# Codexサンドボックスでのrename制限回避のため、書き込み可能な固定パスを使用する
uv_cache_dir := "/home/grace/.codex/uv-cache"
uv_tmp_dir := "/home/grace/.codex/uv-tmp"
web_port := "3001"
api_port := "8000"

@default:
  just --list

setup:
  #!/usr/bin/env bash
  set -euo pipefail
  ./scripts/setup-dev.sh

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
  source ../../scripts/just_env.sh
  load_just_env
  if [[ -z "${UV_CACHE_DIR:-}" ]]; then
    if command -v python3 >/dev/null 2>&1; then
      cache_dir="$(python3 ../../scripts/uv_pick_cache.py "{{uv_cache_dir}}" "$HOME/.cache/uv" "/tmp/uv-cache")" || {
        echo "UV cache でクロスディレクトリrenameが失敗します。Linux ext4上のパスをUV_CACHE_DIRで指定してください。" >&2
        exit 1
      }
    else
      cache_dir="$HOME/.cache/uv"
    fi
  else
    cache_dir="$UV_CACHE_DIR"
  fi
  tmp_dir="${TMPDIR:-$cache_dir/tmp}"
  mkdir -p "$cache_dir" "$tmp_dir"
  TMPDIR="$tmp_dir" UV_CACHE_DIR="$cache_dir" uv sync
  TMPDIR="$tmp_dir" UV_CACHE_DIR="$cache_dir" uv run uvicorn main:app --reload --port "$port"

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
  source ../../scripts/just_env.sh
  load_just_env
  if [[ -z "${UV_CACHE_DIR:-}" ]]; then
    if command -v python3 >/dev/null 2>&1; then
      cache_dir="$(python3 ../../scripts/uv_pick_cache.py "{{uv_cache_dir}}" "$HOME/.cache/uv" "/tmp/uv-cache")" || {
        echo "UV cache でクロスディレクトリrenameが失敗します。Linux ext4上のパスをUV_CACHE_DIRで指定してください。" >&2
        exit 1
      }
    else
      cache_dir="$HOME/.cache/uv"
    fi
  else
    cache_dir="$UV_CACHE_DIR"
  fi
  tmp_dir="${TMPDIR:-$cache_dir/tmp}"
  mkdir -p "$cache_dir" "$tmp_dir"
  TMPDIR="$tmp_dir" UV_CACHE_DIR="$cache_dir" uv sync
  TMPDIR="$tmp_dir" UV_CACHE_DIR="$cache_dir" uv run uvicorn main:app --reload --port "$api_port" &
  popd >/dev/null
  pnpm -C apps/web dev --port "$web_port" &
  wait

zpe-worker:
  #!/usr/bin/env bash
  set -euo pipefail
  ./scripts/run-zpe-worker.sh

zpe-http-worker:
  #!/usr/bin/env bash
  set -euo pipefail
  ./scripts/run-zpe-http-worker.sh

api-test:
  #!/usr/bin/env bash
  set -euo pipefail
  pushd apps/api >/dev/null
  source ../../scripts/just_env.sh
  load_just_env
  if [[ -z "${UV_CACHE_DIR:-}" ]]; then
    if command -v python3 >/dev/null 2>&1; then
      cache_dir="$(python3 ../../scripts/uv_pick_cache.py "{{uv_cache_dir}}" "$HOME/.cache/uv" "/tmp/uv-cache")" || {
        echo "UV cache でクロスディレクトリrenameが失敗します。Linux ext4上のパスをUV_CACHE_DIRで指定してください。" >&2
        exit 1
      }
    else
      cache_dir="$HOME/.cache/uv"
    fi
  else
    cache_dir="$UV_CACHE_DIR"
  fi
  tmp_dir="${TMPDIR:-$cache_dir/tmp}"
  mkdir -p "$cache_dir" "$tmp_dir"
  link_mode="${UV_LINK_MODE:-copy}"
  uv_flags=()
  if [[ "${UV_NO_CACHE:-}" == "1" ]]; then
    uv_flags+=(--no-cache)
  fi
  TMPDIR="$tmp_dir" UV_CACHE_DIR="$cache_dir" UV_LINK_MODE="$link_mode" uv sync "${uv_flags[@]}"
  TMPDIR="$tmp_dir" UV_CACHE_DIR="$cache_dir" UV_LINK_MODE="$link_mode" uv run "${uv_flags[@]}" pytest
  popd >/dev/null

api-ruff:
  #!/usr/bin/env bash
  set -euo pipefail
  pushd apps/api >/dev/null
  source ../../scripts/just_env.sh
  load_just_env
  if [[ -z "${UV_CACHE_DIR:-}" ]]; then
    if command -v python3 >/dev/null 2>&1; then
      cache_dir="$(python3 ../../scripts/uv_pick_cache.py "{{uv_cache_dir}}" "$HOME/.cache/uv" "/tmp/uv-cache")" || {
        echo "UV cache でクロスディレクトリrenameが失敗します。Linux ext4上のパスをUV_CACHE_DIRで指定してください。" >&2
        exit 1
      }
    else
      cache_dir="$HOME/.cache/uv"
    fi
  else
    cache_dir="$UV_CACHE_DIR"
  fi
  tmp_dir="${TMPDIR:-$cache_dir/tmp}"
  mkdir -p "$cache_dir" "$tmp_dir"
  link_mode="${UV_LINK_MODE:-copy}"
  uv_flags=()
  if [[ "${UV_NO_CACHE:-}" == "1" ]]; then
    uv_flags+=(--no-cache)
  fi
  TMPDIR="$tmp_dir" UV_CACHE_DIR="$cache_dir" UV_LINK_MODE="$link_mode" uv sync "${uv_flags[@]}"
  TMPDIR="$tmp_dir" UV_CACHE_DIR="$cache_dir" UV_LINK_MODE="$link_mode" uv run "${uv_flags[@]}" ruff check .
  popd >/dev/null

api-mypy:
  #!/usr/bin/env bash
  set -euo pipefail
  pushd apps/api >/dev/null
  source ../../scripts/just_env.sh
  load_just_env
  if [[ -z "${UV_CACHE_DIR:-}" ]]; then
    if command -v python3 >/dev/null 2>&1; then
      cache_dir="$(python3 ../../scripts/uv_pick_cache.py "{{uv_cache_dir}}" "$HOME/.cache/uv" "/tmp/uv-cache")" || {
        echo "UV cache でクロスディレクトリrenameが失敗します。Linux ext4上のパスをUV_CACHE_DIRで指定してください。" >&2
        exit 1
      }
    else
      cache_dir="$HOME/.cache/uv"
    fi
  else
    cache_dir="$UV_CACHE_DIR"
  fi
  tmp_dir="${TMPDIR:-$cache_dir/tmp}"
  mkdir -p "$cache_dir" "$tmp_dir"
  link_mode="${UV_LINK_MODE:-copy}"
  uv_flags=()
  if [[ "${UV_NO_CACHE:-}" == "1" ]]; then
    uv_flags+=(--no-cache)
  fi
  TMPDIR="$tmp_dir" UV_CACHE_DIR="$cache_dir" UV_LINK_MODE="$link_mode" uv sync "${uv_flags[@]}"
  TMPDIR="$tmp_dir" UV_CACHE_DIR="$cache_dir" UV_LINK_MODE="$link_mode" uv run "${uv_flags[@]}" mypy .
  popd >/dev/null
