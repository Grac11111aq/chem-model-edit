# Tasks: ZPE 計算サーバ分離

**Input**: `specs/003-zpe-remote-compute/spec.md` / `specs/003-zpe-remote-compute/plan.md`
**Tests**: pytest（バックエンド/ストアの単体テストを優先、E2E は環境確立後に追加）

**Format**: `[ID] [P?] [Story] Description`
**Stories**: US1=委譲, US2=結果共有, US3=モック

---

## Phase 1: Setup

- [ ] T200 [P] [US1] ZPE 設定に `compute_mode` / `result_store` / `auth` を追加
- [ ] T201 [P] [US2] `ResultStore` インターフェースと Redis 実装を作成

## Phase 2: Core Implementation

- [ ] T210 [P] [US1] `remote-queue` バックエンド（RQ への enqueue + meta 更新）
- [ ] T211 [P] [US2] 結果取得 API が Redis ストア経由で返すように変更
- [ ] T212 [US3] `mock` バックエンド（決定論的ダミー結果）
- [ ] T213 [US1/2] エラーハンドリング（計算サーバー/ストア障害の明確化）
- [ ] T214 [US1] 計算サーバ登録トークン（MVP: 管理者発行、将来拡張前提）

## Phase 3: Tests & Validation

- [ ] T220 [P] [US1/2] Redis ストアの単体テスト
- [ ] T221 [US3] モックモードのAPIテスト
- [ ] T222 [P] [US1/2/3] pytest / mypy / ruff 実行

## Phase 4: Docs & Ops

- [x] T230 [P] [US1/2] compute-plane セットアップガイドと env テンプレート整備
- [x] T231 [P] [US1] worker 起動スクリプト（+ just レシピ）追加
- [x] T232 [P] [US1/2] README / setup ガイドの導線追加
