#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

echo "==> Node + pnpm"
if ! command -v corepack >/dev/null 2>&1; then
  echo "corepack not found. Install Node.js (with Corepack) first." >&2
  exit 1
fi
corepack enable
corepack prepare pnpm@10.27.0 --activate
pnpm install

if [[ "${SETUP_SKIP_APPROVE:-}" != "1" ]]; then
  if grep -q "@parcel/watcher" pnpm-workspace.yaml; then
    echo "pnpm approve-builds already applied for @parcel/watcher."
  else
    echo "==> Approve build scripts (select @parcel/watcher)"
    pnpm approve-builds
  fi
fi

if command -v uv >/dev/null 2>&1; then
  echo "==> Python + API deps (uv)"
  uv python install 3.13
  (cd apps/api && uv sync)
else
  echo "uv not found. Install uv to set up Python deps." >&2
fi

echo "Setup done."
