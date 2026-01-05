# Tasks: 計算化学向け構造編集Webアプリ

**Input**: Design documents from `specs/001-chem-model-webapp/`
**Prerequisites**: `specs/001-chem-model-webapp/plan.md`, `specs/001-chem-model-webapp/spec.md`, `specs/001-chem-model-webapp/ui-ux.md`

**Tests**: pytest を導入済みのため、APIテストを優先して段階的に追加する。

**Format**: `[ID] [P?] [Story] Description`  
**Stories**: US1=読み込み/編集, US2=比較/移植, US3=スーパーセル/共有

---

## Phase 1: Setup (Shared Infrastructure)

- [x] T001 [P] [US1] Python API プロジェクト初期化（`apps/api/` + uv）
- [x] T002 [P] [US1] `ruff`/`mypy`/`pytest` を `apps/api/pyproject.toml` に追加
- [x] T003 [P] [US1] pnpm workspace を用意し `pnpm-workspace.yaml` を追加
- [x] T004 [P] [US1] `apps/web` を TanStack Start で初期化
- [x] T005 [P] [US1] フロント共通の lint/format 設定を `apps/web/` に反映
- [x] T006 [P] [US1] `.gitignore` に `.venv/` など生成物を追加

---

## Phase 2: Foundational (Blocking)

- [x] T010 [US1] FastAPI アプリ骨組み作成（`apps/api/main.py`）
- [x] T011 [P] [US1] Pydantic モデル定義（`apps/api/models.py`）
- [x] T012 [US1] .in パースサービス（ASE/pymatgen）実装（`apps/api/services/parse.py`）
- [x] T013 [US1] エクスポートサービス（.in書き戻し）実装（`apps/api/services/export.py`）
- [x] T014 [US1] API ルーティング（`/parse`, `/export`）を実装（`apps/api/main.py`）
- [x] T015 [US1] フロント基盤: ルーティングとレイアウトを実装（`apps/web/src/routes/*`, `apps/web/src/components/layout/*`）
- [x] T016 [US1] Mol* Viewer コンポーネント基盤（`apps/web/src/components/molstar/*`）

---

## Phase 3: User Story 1 - .in の読み込み/編集 (P1) 🎯

**Goal**: .in から構造を読み込み、テーブル編集と3Dが同期する

- [x] T020 [P] [US1] API: 解析結果の型定義を共通化（`packages/shared/src/types.ts`）
- [x] T021 [US1] フロント: インポートUI（`apps/web/src/routes/editor.tsx`）
- [x] T022 [US1] フロント: Atom Table（編集可能）実装（`apps/web/src/routes/editor.tsx`）
- [x] T023 [US1] フロント: 3D表示に座標編集を反映（`apps/web/src/routes/editor.tsx`, `apps/web/src/components/molstar/*`）
- [x] T024 [US1] エクスポートUI（.in 保存）実装（`apps/web/src/routes/editor.tsx`）

### Tests (US1)
- [x] T025 [P] [US1] API: `/parse` の単体テスト（`apps/api/tests/test_parse.py`）
- [x] T026 [P] [US1] API: `/export` の単体テスト（`apps/api/tests/test_export.py`）

---

## Phase 4: User Story 2 - 比較・部分移植 (P2)

**Goal**: 複数構造を並べて比較/移植/整列できる

- [x] T030 [US2] フロント: 構造リストUI（任意数管理）（`apps/web/src/routes/editor.tsx`）
- [x] T031 [US2] フロント: 選択/コピー/貼り付け操作（`apps/web/src/routes/editor.tsx`）
- [x] T032 [US2] フロント: 距離テーブル（`apps/web/src/routes/editor.tsx`）
- [x] T033 [US2] フロント: 整列/シフト操作（`apps/web/src/components/compare/align.ts`）

---

## Phase 5: User Story 3 - スーパーセル/共有 (P3)

**Goal**: スーパーセル生成と共有（単一HTML）

- [x] T040 [US3] API: スーパーセル生成サービス（`apps/api/services/supercell.py`）
- [x] T041 [US3] フロント: スーパーセル画面UI（`apps/web/src/routes/supercell.tsx`）
- [x] T042 [US3] フロント: 単一HTMLエクスポート（`apps/web/src/components/share/html-export.tsx`）
- [x] T043 [US3] 共有リンク方式の設計メモを追加（`specs/001-chem-model-webapp/share-notes.md`）

---

## Phase 6: Polish & Cross-Cutting

- [x] T050 [P] UI/UX の細部調整（ショートカット、アクセシビリティ）
- [x] T051 [P] パフォーマンス最適化（大規模構造の表示）
- [x] T052 [P] ドキュメント更新（`specs/001-chem-model-webapp/quickstart.md`）
