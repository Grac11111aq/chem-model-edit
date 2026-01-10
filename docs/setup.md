# Setup (local development)

This guide covers the local dev environment for `chem-model-edit`.
Quantum ESPRESSO (QE) runtime is **not** required here.

## Prerequisites
- Git
- Node.js with Corepack enabled
- pnpm 10.27.0
- Python >= 3.13 (use `uv`)
- Optional: `just` (recommended for `just dev`)

If you want `just`:
- With Rust/Cargo installed: `cargo install just`
- Or via your OS package manager

## Quick start (script)
```bash
./scripts/setup-dev.sh
```
Notes:
- The script runs `pnpm approve-builds` and asks you to approve `@parcel/watcher`.
- `pnpm approve-builds` will update `pnpm-workspace.yaml`.
- Set `SETUP_SKIP_APPROVE=1` to skip the approval step.

## Manual setup
```bash
# Node + pnpm
corepack enable
corepack prepare pnpm@10.27.0 --activate
pnpm install

# Approve build scripts (select @parcel/watcher)
pnpm approve-builds

# Python + API deps
uv python install 3.13
cd apps/api
uv sync
```

## Run
```bash
# Web
pnpm -C apps/web dev --port 3000

# API
cd apps/api
uv run uvicorn main:app --reload --port 8000
```

Run both together (requires `just`):
```bash
just dev
```
If `just` complains about `XDG_RUNTIME_DIR`, use the wrapper:
```bash
./scripts/just dev
```

## Quality checks
```bash
pnpm -C apps/web typecheck
pnpm -C apps/web lint
pnpm -C apps/web test

cd apps/api
uv run pytest
uv run mypy .
```

## Troubleshooting
- In some sandboxed environments, `uv sync` can fail with `EXDEV` (cross-device rename).
  Use a non-sandboxed environment or set `UV_CACHE_DIR`/`TMPDIR` to a path on the same
  filesystem. The `Justfile` already pins these paths for the Codex sandbox.
