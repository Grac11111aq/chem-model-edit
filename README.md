# chem-model-edit

English | [Japanese](./README.ja.md)

Note: the Japanese README is a short summary; this English version is the canonical reference.

Web app for editing and visualizing chemical structures (Quantum ESPRESSO `.in`).

## Project goal
Build a web-based workflow to view, edit, compare, and share QE `.in` structures.
For the detailed scope and acceptance criteria, see `specs/001-chem-model-webapp/spec.md`.

## Scope (per spec)
- Parse QE `.in` files (starting with species + coordinates; lattice support is incremental).
- 3D visualization with Mol* (ball-and-stick by default).
- Table-based editing that stays in sync with 3D view.
- Multi-structure workflows such as transfer/compare (planned).
- Export to QE `.in` and shareable single HTML (planned).

## Repository layout
- `apps/web`: TanStack Start web app (frontend).
- `apps/api`: FastAPI backend.
- `packages/shared`: Shared types.
- `specs/`: Spec/plan/task documents.
- `samples/qe-in`: Sample QE `.in` inputs.
- `ref-legacy`: Legacy reference implementation.

## Prerequisites
- Git
- Node.js (Corepack) + pnpm 10.27.0
- Python >= 3.13 + uv
- Optional: just

## Setup
```bash
git clone git@github.com:Grac11111aq/chem-model-edit.git
cd chem-model-edit

corepack enable
corepack prepare pnpm@10.27.0 --activate
pnpm install
```

## Development
All commands below run from the repo root unless noted.

### Web
```bash
pnpm dev
# or
pnpm -C apps/web dev
```
Vite uses port 3000 by default (unless overridden).

### API
```bash
cd apps/api
uv sync
uv run uvicorn main:app --reload --port 8000
```

### Web + API together
If you have `just` installed, it starts both and auto-selects free ports:
```bash
just dev
```
You can override ports:
```bash
WEB_PORT=3001 API_PORT=8000 just dev
```

## Quality checks
### Web
```bash
pnpm typecheck
pnpm lint
pnpm test
```

### API
```bash
cd apps/api
uv run pytest
uv run mypy .
```

## Specs and samples
- Spec/Plan/Tasks: `specs/001-chem-model-webapp/`
- Sample inputs: `samples/qe-in/`
