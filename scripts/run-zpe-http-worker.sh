#!/usr/bin/env bash
set -euo pipefail

root_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

env_file="${ZPE_ENV_FILE:-}"
if [[ -z "$env_file" && -f "$root_dir/apps/api/.env.compute" ]]; then
  env_file="$root_dir/apps/api/.env.compute"
fi

if [[ -n "$env_file" ]]; then
  if [[ ! -f "$env_file" ]]; then
    echo "Env file not found: $env_file" >&2
    exit 1
  fi
  set -a
  # shellcheck disable=SC1090
  source "$env_file"
  set +a
  export ZPE_ENV_FILE="$env_file"
fi

cd "$root_dir/apps/api"

if [[ "${ZPE_WORKER_SYNC:-}" == "1" ]]; then
  uv sync
fi

exec uv run python -m services.zpe.http_worker "$@"
