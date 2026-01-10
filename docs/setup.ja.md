# セットアップ（ローカル開発）

`chem-model-edit` のローカル開発環境を整える手順です。
Quantum ESPRESSO (QE) の実行環境は **不要** です。

## 前提
- Git
- Node.js（Corepack 有効化）
- pnpm 10.27.0
- Python >= 3.13（`uv` を推奨）
- 任意: `just`（`just dev` を使う場合）

`just` を入れる場合:
- Rust/Cargo があるなら: `cargo install just`
- もしくは OS のパッケージマネージャ

## クイックスタート（スクリプト）
```bash
./scripts/setup-dev.sh
```
補足:
- `pnpm approve-builds` が走るので、`@parcel/watcher` を承認してください。
- `pnpm approve-builds` により `pnpm-workspace.yaml` が更新されます。
- `SETUP_SKIP_APPROVE=1` を指定すると承認手順をスキップできます。

## 手動セットアップ
```bash
# Node + pnpm
corepack enable
corepack prepare pnpm@10.27.0 --activate
pnpm install

# build scripts を承認（@parcel/watcher を選択）
pnpm approve-builds

# Python + API 依存関係
uv python install 3.13
cd apps/api
uv sync
```

## 起動
```bash
# Web
pnpm -C apps/web dev --port 3000

# API
cd apps/api
uv run uvicorn main:app --reload --port 8000
```

同時起動（`just` が必要）:
```bash
just dev
```
`XDG_RUNTIME_DIR` のエラーが出る場合はラッパーを使用:
```bash
./scripts/just dev
```

## 品質チェック
```bash
pnpm -C apps/web typecheck
pnpm -C apps/web lint
pnpm -C apps/web test

cd apps/api
uv run pytest
uv run mypy .
```

## トラブルシュート
- サンドボックス環境では `uv sync` が `EXDEV`（cross-device rename）で失敗する場合があります。
  同一ファイルシステム上の `UV_CACHE_DIR` / `TMPDIR` を指定するか、非サンドボックス環境で実行してください。
  Codex 向けには `Justfile` が固定パスを使うようにしてあります。
