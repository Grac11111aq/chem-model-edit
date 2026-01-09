# Tasks: ZPEGUI バックエンド移植

**Input**: `specs/002-zpegui-backend/spec.md` / `specs/002-zpegui-backend/plan.md`
**Tests**: pytest を優先（解析ユーティリティとAPIの単体テスト）

**Format**: `[ID] [P?] [Story] Description`
**Stories**: US1=解析, US2=ジョブ投入, US3=状態/結果取得

---

## Phase 1: Setup

- [x] T100 [P] [US2] API 依存追加（`rq`, `redis`, `pydantic-settings` 等）
- [x] T101 [P] [US2] ZPE 設定管理（Redis URL / 作業ディレクトリ / QE 実行系）
- [x] T102 [US2] RQ worker エントリポイント作成

## Phase 2: Core Implementation

- [x] T110 [P] [US1] QE .in の固定原子抽出・解析ユーティリティ作成
- [x] T111 [P] [US1] `/calc/zpe/parse` の request/response モデル追加
- [x] T112 [P] [US2] ZPE ジョブ投入 API（`/calc/zpe/jobs`）
- [x] T113 [US2] RQ ジョブ本体（Vibrations / キャッシュ健全性 / environ 配布 / ZPE算出）
- [x] T114 [US3] ジョブ状態取得 API（`/calc/zpe/jobs/{id}`）
- [x] T115 [US3] 結果取得/CSV DL API（`/calc/zpe/jobs/{id}/result` / `/calc/zpe/jobs/{id}/files`）

## Phase 3: Tests & Validation

- [x] T120 [P] [US1] 解析ユーティリティの単体テスト
- [x] T121 [US3] ジョブメタ/状態レスポンスの単体テスト（fakeredis を検討）
- [ ] T122 [P] [US1/2/3] pytest / mypy / ruff 実行
