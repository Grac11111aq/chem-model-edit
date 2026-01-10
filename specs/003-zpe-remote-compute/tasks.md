# Tasks: ZPE 計算サーバ分離

**Input**: `specs/003-zpe-remote-compute/spec.md` / `specs/003-zpe-remote-compute/plan.md`
**Tests**: pytest（バックエンド/ストアの単体テストを優先）

**Format**: `[ID] [P?] [Story] Description`
**Stories**: US1=委譲, US2=結果共有, US3=モック

---

## Phase 1: Setup

- [x] T200 [P] [US1] ZPE 設定に `compute_mode` / `result_store` / `auth` を追加
- [x] T201 [P] [US2] `ResultStore` インターフェースと Redis 実装を作成

## Phase 2: Core Implementation

- [x] T210 [P] [US1] `remote-queue` バックエンド（RQ への enqueue + meta 更新）
- [x] T211 [P] [US2] 結果取得 API が Redis ストア経由で返すように変更
- [x] T212 [US3] `mock` バックエンド（決定論的ダミー結果）
- [x] T213 [US1/2] エラーハンドリング（計算サーバー/ストア障害の明確化）
- [x] T214 [US1] 計算サーバ登録トークン（MVP: 管理者発行、将来拡張前提）

## Phase 3: Tests & Validation

- [x] T220 [P] [US1/2] Redis ストアの単体テスト
- [x] T221 [US3] モックモードのAPIテスト
- [x] T222 [P] [US1/2/3] pytest / mypy / ruff 実行
