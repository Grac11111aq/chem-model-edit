#!/usr/bin/env bash
set -euo pipefail

load_just_env() {
  local env_files=()
  if [[ -n "${JUST_ENV_FILE:-}" ]]; then
    env_files+=("$JUST_ENV_FILE")
  fi
  env_files+=(".just.env" "${XDG_CONFIG_HOME:-$HOME/.config}/chem-model-edit/just.env")

  for env_file in "${env_files[@]}"; do
    if [[ -f "$env_file" ]]; then
      set -a
      # shellcheck disable=SC1090
      source "$env_file"
      set +a
    fi
  done

  if [[ -z "${TMPDIR:-}" && -n "${UV_TMP_DIR:-}" ]]; then
    TMPDIR="$UV_TMP_DIR"
    export TMPDIR
  fi
}
