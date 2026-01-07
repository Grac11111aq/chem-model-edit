# Feature Specification: README改善

**Feature Branch**: `[docs/readme-refresh]`  
**Created**: 2026-01-07  
**Status**: Draft  
**Input**: ユーザー要望「Readmeの改善」

## 背景 / 目的
- 現状のREADMEはセットアップ/起動の最低限情報が中心で、目的・スコープ・構成が把握しづらい。
- `apps/web` のREADMEはテンプレート文面で、トップREADMEにプロジェクト固有の導線を整理したい。

## スコープ
- READMEの構成整理と情報の拡充（概要/スコープ/構成/開発/品質/参照リンク）。
- 既存コマンドの正確性チェック（`package.json` / `Justfile` / `apps/api/pyproject.toml` に合わせる）。
- `README.ja.md` を追加し、日本語版は無理な直訳を避けた要点版として整備する。
- README内に言語切替リンクを設ける。

## 非ゴール
- 実装機能の追加・変更。
- CI/CDの導入（別Issueで管理）。
- アプリUIの多言語化。
- 全面的なドキュメント再構成。

## 成功条件
- READMEだけで「何を作っているか / どこを見れば良いか / どう動かすか」が3分以内に把握できる。
- 記載コマンドが実際の設定と矛盾しない。

## 受け入れ基準
- READMEに以下のセクションが含まれる：概要、スコープ/機能、リポジトリ構成、開発手順（Web/API/同時起動）、品質チェック、参考リンク。
- 主要コマンドとして `pnpm dev`, `pnpm -C apps/web dev`, `uv run uvicorn ...`, `just dev` が記載される。
- `specs/` と `samples/qe-in` への導線がある。
- `README.md` と `README.ja.md` が相互リンクされる。
- 日本語版は無理な直訳を避け、必要に応じて英語用語を併記する。

## 検証方法
- README記述と `package.json`, `Justfile`, `apps/api/pyproject.toml` の整合性を確認する。
