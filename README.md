# chem-model-edit

Web app for editing and visualizing chemical structures (Quantum ESPRESSO `.in`).

## Requirements
- Git
- Node.js (Corepack) + pnpm 10.27.0
- Python >= 3.13 + uv
- Optional: just

## Setup
```bash
git clone https://github.com/Grac11111aq/chem-model-edit.git
cd chem-model-edit

corepack enable
corepack prepare pnpm@10.27.0 --activate
pnpm install
```

## Run (web)
```bash
pnpm dev
# or
pnpm -C apps/web dev
```
Default Vite port is 3000.

## Run (api)
```bash
cd apps/api
uv sync
uv run uvicorn main:app --reload --port 8000
```

## Run both (optional)
If you have `just` installed, it starts both processes and finds free ports:
```bash
just dev
```
You can override ports:
```bash
WEB_PORT=3001 API_PORT=8000 just dev
```

## Quality checks
```bash
pnpm typecheck
pnpm lint
pnpm test
```
