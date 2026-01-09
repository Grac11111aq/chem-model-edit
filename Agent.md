# Codex Agent Context (chem-model-edit)

## 目的 / 概要
計算化学向けの構造可視化・編集ツールをWebアプリとして提供する。複数構造の並列表示、部分的な構造移植、比較・整列、エクスポート共有を支援する。

## スコープ（現時点）
- Quantum ESPRESSOの`.in`を対象に構造を読み取り・編集・書き出し
- 初期は座標と原子種のみを扱う
- 共有は「単一HTML」(NGL埋め込み) を提供。将来的にサーバ保存リンクを優先
- 3D表示は Mol* を使い、**デフォルトはボール＆スティック**

## 技術スタック
- Frontend: TanStack Start + shadcn/ui + Mol*
- Backend: FastAPI + ASE / pymatgen
- JS/TS: pnpm
- Python: uv (+ ruff, mypy, pytest)

## リポジトリ構成
- `apps/web`: Webアプリ (TanStack Start)
- `apps/api`: FastAPI API
- `packages/shared`: 共通型
- `specs/001-chem-model-webapp`: Spec/Plan/Task
- `ref-legacy`: 旧実装 (参照用)
- `samples/qe-in`: テスト用の Quantum ESPRESSO `.in` サンプル
- `uploads`: ユーザーが外部ドキュメント/データ/サンプルをアップロードする場所

## 重要なUI/機能 (web)
- 複数構造管理、選択コピー/貼り付け、距離表示、シフト/整列
- Import: ファイル/クリップボード
- Export: QE `.in` / Share HTML
- Supercell生成: API呼び出し + Mol*プレビュー

## エージェント作業方針
- **日本語で報告**
- Spec Kit ワークフロー順守:
  - Spec → Plan → Task を先に整備
  - 変更が想定外ならSpecを更新して再承認
- 既存の仕様/タスクに従う
- `uploads` 配下のユーザーアップロード済みデータは、作業上必要であれば移動・整理・削除してよい

## 開発/検証コマンド
- Web:
  - `pnpm -C apps/web dev --port 3001`
  - `XDG_RUNTIME_DIR=.xdg-runtime just web` (justが権限エラーになる場合)
  - `pnpm -C apps/web typecheck`
- API:
  - `uv run uvicorn apps.api.main:app --reload --port 8000`
  - `uv run pytest`
  - `uv run mypy .`

## 既知の注意点
- Mol*は`Viewer.create`を使い、PDB読み込み後に`ball-and-stick`表現を追加
- Import/Exportの入出力は`.in` (QE) と座標/原子種に限定
- `apps/web/src/components/molstar/MolstarViewer.tsx` が表示の中心
- pnpmストアは `/home/grace/projects/chem-model-edit/.pnpm-store` を共通利用する方針（`.pnpm-store` はgit ignore）
- Codexサンドボックス（workspace-write）では `uv sync` が rename 制限で失敗するため、uv 実行時は昇格（on-request）か danger-full-access が必要

## 進め方のベストプラクティス
- 変更前に影響範囲と対象ファイルを明確化
- 破壊的変更は避け、既存実装との互換性を維持
- 必要に応じてテスト/型チェックを実行
- 実装結果はSpec/Plan/Taskへ反映
- Chrome-devtools-mcpを使ってブラウザでの動作検証
- 作業の区切りのたびにgitを利用してください
- 編集後は型チェック，Lint, テストを行って
