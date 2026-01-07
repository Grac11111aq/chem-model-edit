# Issue: CI/CD セットアップ

**Created**: 2026-01-07  
**Status**: Open  
**Owner**: TBD

## 背景
- 現状はCI/CD未導入で、PR時の品質チェックが手動。
- Web/APIの両方に検証コマンドがあり、自動化の優先度が高い。

## 目的
- PR/Pushで最低限の品質チェックを自動化する。
- 失敗時に明確なステータスが出るようにする。

## 想定スコープ（案）
- GitHub ActionsでPR/Pushトリガーを設定。
- Web: `pnpm -C apps/web typecheck`, `pnpm -C apps/web lint`, `pnpm -C apps/web test`。
- API: `uv sync`, `uv run pytest`, `uv run mypy .`（必要なら ruff を追加）。
- 依存キャッシュ（pnpm/uv）を導入。

## 非ゴール
- デプロイやリリース自動化。
- 本番環境の運用設定。

## 受け入れ条件
- PR作成時にCIが自動実行される。
- 失敗時にGitHub上で赤ステータスになる。
- 主要チェックが再現可能なコマンドとして文書化される。

## 未決事項 / Open Questions
- Node/Pythonのバージョン固定（例: Node 20+ / Python 3.13）。
- テスト時間と並列化（Web/APIを分けるか）。
- ruff/mypyの設定ファイル追加が必要か。
- ActionsバッジをREADMEに追加するタイミング。

## 次アクション
- ステークホルダーの合意を得てSpec/Plan化する。
