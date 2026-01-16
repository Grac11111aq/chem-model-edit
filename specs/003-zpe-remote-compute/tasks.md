# Tasks: ZPE 計算サーバ分離

**Input**: `specs/003-zpe-remote-compute/spec.md` / `specs/003-zpe-remote-compute/plan.md`
**Tests**: pytest（バックエンド/ストアの単体テストを優先、E2E は環境確立後に追加）

**Format**: `[ID] [P?] [Story] Description`
**Stories**: US1=委譲, US2=結果共有, US3=モック, US4=HTTP仲介

---

## Phase 1: Setup

- [x] T200 [P] [US1] ZPE 設定に `compute_mode` / `result_store` / `auth` を追加
- [x] T201 [P] [US2] `ResultStore` インターフェースと Redis 実装を作成

## Phase 2: Core Implementation (HTTP仲介)

- [x] T210 [P] [US4] `remote-http` enqueue（payload保存 + queue投入）
- [x] T211 [P] [US4] worker 登録トークン + worker_token 発行/失効
- [x] T212 [P] [US4] lease 取得 API + 期限切れ再投入
- [x] T213 [P] [US4] result/failed API + retry/backoff/DLQ
- [ ] T214 [P] [US4] HTTP ワーカーポーリング実装
- [ ] T215 [US4] worker 起動・設定フロー（CLI/スクリプト）

## Phase 3: Tests & Validation

- [x] T220 [P] [US4] HTTP enqueue/lease/result/failed の単体テスト
- [ ] T221 [P] [US4] HTTP ワーカーの結合テスト
- [ ] T222 [US3] モックモードのAPIテスト
- [ ] T223 [P] [US1/2/3/4] pytest / mypy / ruff 実行

## Phase 4: Docs & Ops

- [ ] T230 [P] [US4] HTTP ワーカーセットアップガイド（Upstash非公開前提）
- [ ] T231 [P] [US4] H2最小サンプル入力と送信スクリプト追加
- [ ] T232 [P] [US4] E2E 検証手順をSpec/Planに反映

### Notes
- T221 は T214 完了後に着手（HTTP ワーカー実装が前提）
- T223 は control-plane 実装の早期検証として先行実行可能
- T231 はワーカー実装とは独立なので先に着手可能
