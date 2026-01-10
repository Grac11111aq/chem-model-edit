# chem-model-edit

[English](./README.md) | 日本語

Quantum ESPRESSO `.in` の構造を可視化・編集するWebアプリ。  
※日本語版は要点のみです。詳細/最新は英語版を参照してください。

## 目的 / ゴール
- QE `.in` 構造の view / edit / compare / share をブラウザで完結させる。
- 詳細なスコープと受け入れ条件は `specs/001-chem-model-webapp/spec.md` を参照。

## スコープ（Specに基づく）
- QE `.in` の読み込み（まずは species + coordinates から）。
- Mol* による 3D visualization（デフォルトは ball-and-stick）。
- テーブル編集と 3D 表示の同期。
- 複数構造の transfer / compare（planned）。
- QE `.in` 出力と shareable single HTML（planned）。

## リポジトリ構成
- `apps/web`: TanStack Start web app（frontend）。
- `apps/api`: FastAPI backend。
- `packages/shared`: Shared types。
- `specs/`: Spec/plan/task。
- `samples/qe-in`: QE `.in` サンプル。
- `ref-legacy`: 旧実装（参照用）。

## 前提
- Git
- Node.js (Corepack) + pnpm 10.27.0
- Python >= 3.13 + uv
- Optional: just

## セットアップ
詳細は `docs/setup.ja.md` を参照。
```bash
git clone git@github.com:Grac11111aq/chem-model-edit.git
cd chem-model-edit

corepack enable
corepack prepare pnpm@10.27.0 --activate
pnpm install
# @parcel/watcher を承認
pnpm approve-builds
```

## 開発
※以下のコマンドは、特記がない限り repo root から実行します。

### Web
```bash
pnpm dev
# or
pnpm -C apps/web dev
```
Vite のデフォルトは 3000（上書きがない場合）。

### API
```bash
cd apps/api
uv sync
uv run uvicorn main:app --reload --port 8000
```

### Web + API 同時起動
`just` がある場合は両方を起動し、空きポートを自動選択します。
```bash
just dev
```
ポート指定:
```bash
WEB_PORT=3001 API_PORT=8000 just dev
```

## 品質チェック
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

## 参考
- Spec/Plan/Tasks: `specs/001-chem-model-webapp/`
- サンプル入力: `samples/qe-in/`
