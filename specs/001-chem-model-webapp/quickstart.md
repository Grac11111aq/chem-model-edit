# Quickstart: Chem Model Edit

## 前提
- Node.js + pnpm
- Python 3.13+（uv で環境管理）

## セットアップ

### Frontend
```bash
pnpm install
pnpm -C apps/web dev
```

### Backend
```bash
cd apps/api
uv sync
uv run uvicorn main:app --reload
```

## 動作確認

### 1) Editor
1. `http://localhost:3000/editor` を開く
2. `.in` を貼り付けて Parse
3. Atom Table で編集 → Mol* に反映
4. Export .in / Share HTML で出力

### 2) Supercell
1. `http://localhost:3000/supercell` を開く
2. A/B のXYZとシーケンスを入力
3. Generate → Mol* プレビュー

## API エンドポイント
- `POST /parse` `.in` → atoms
- `POST /export` atoms → `.in`（ATOMIC_POSITIONSのみ）
- `POST /supercell` A/B + sequence + lattice → atoms

## 注意
- Share HTML は NGL のCDNを使用（オフラインでは表示不可）
- 大規模構造は描画に時間がかかるため、編集後に反映まで遅延が発生する場合あり
