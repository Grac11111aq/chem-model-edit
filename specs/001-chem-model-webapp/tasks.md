# Tasks: 計算化学向け構造編集Webアプリ

**Input**: Design documents from `specs/001-chem-model-webapp/`
**Prerequisites**: `specs/001-chem-model-webapp/plan.md`, `specs/001-chem-model-webapp/spec.md`, `specs/001-chem-model-webapp/ui-ux.md`

**Tests**: pytest を導入済みのため、APIテストを優先して段階的に追加する。

**Format**: `[ID] [P?] [Story] Description`  
**Stories**: US1=読み込み/編集, US2=表面転写/距離QC, US3=複合スーパーセル, US4=格子変換, US5=共有/重ね表示

**Note**: 共有は当面単一HTMLで実装し、重ね表示/可視/透明度を再現する。

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
- [x] T027 [US1] 出力設定画面（入力形式/単位保持がデフォルト、将来拡張の枠を用意）

### Tests (US1)
- [x] T025 [P] [US1] API: `/parse` の単体テスト（`apps/api/tests/test_parse.py`）
- [x] T026 [P] [US1] API: `/export` の単体テスト（`apps/api/tests/test_export.py`）

---

## Phase 4: User Story 2 - 表面転写/距離QC (P2)

**Goal**: 2構造(A/B)の表面転写と対応原子距離のQCができる

- [x] T030 [US2] フロント: 比較対象の選択UI（A/Bを中心に切替）
- [x] T031 [US2] フロント: 表面領域の選択UI（手動選択をベース）
- [x] T032 [US2] 表面転写（インデックス対応、Cartesian座標の上書き）
- [x] T033 [US2] 対応原子距離レポート（PBC最小像を適用）
- [x] T034 [US2] 事前検証（原子数非互換の警告）
- [x] T035 [US2] フロント: 距離テーブルUI
- [x] T036 [US5] フロント: 重ね表示ON/OFFと構造ごとの表示/透明度をMol*へ反映（`apps/web/src/components/molstar/MolstarViewer.tsx`, `apps/web/src/routes/editor.tsx`）
- [x] T037 [US5] フロント: 構造ごとの表示/透明度UIを追加（`apps/web/src/routes/editor.tsx`）

---

## Phase 4b: User Story 2b - Δ移植 (P2)

**Goal**: 小スラブ出力の最終座標から変位を算出し、大スラブへ移植する

- [x] T038 [US2b] API: Δ移植サービス実装（入力: small_in/out/large_in, 出力: transplanted .in）
- [x] T039 [US2b] API: ルート追加（`/transplant/delta`）
- [x] T040 [US2b] API: Δ移植のテスト追加（正常系/異常系）
- [x] T041 [US2b] フロント: 専用ページUI（3入力 + 実行 + 結果表示）
- [x] T042 [US2b] フロント: 結果のコピー/ダウンロード

---

## Phase 5: User Story 3 - 複合スーパーセル (P3)

**Goal**: タイルパターンに従って複合スーパーセルを生成する

- [x] T040 [US3] API: 複合スーパーセル生成（タイルパターン入力）
- [x] T041 [US3] フロント: タイルパターン入力UI（2D配列 + GUIの下地）
- [x] T042 [US3] フロント: 生成結果のプレビュー
- [x] T043 [US3] 重複/衝突チェック（許容誤差付き、任意）

---

## Phase 6: User Story 4 - 格子変換 (P3)

**Goal**: 格子定数とセルベクトルを相互変換できる

- [x] T060 [US4] API: 格子変換サービス（a,b,c,α,β,γ ⇄ cell vectors）
- [x] T061 [US4] フロント: 変換UI（単位切替: alat/bohr/angstrom）
- [x] T062 [US4] 変換の単体テスト追加

---

## Phase 7: User Story 5 - 共有/重ね表示 (P3)

**Goal**: 重ね表示の状態を共有HTMLで再現できる

- [x] T070 [US5] 共有HTML: 重ね表示/表示/透明度を再現（`apps/web/src/components/share/html-export.tsx`）
- [x] T071 [US5] 共有HTMLのエラーハンドリング整備（空データ時の通知）

---

## Phase 8: Polish & Cross-Cutting

- [x] T050 [P] UI/UX の細部調整（ショートカット、アクセシビリティ）
- [x] T051 [P] パフォーマンス最適化（大規模構造の表示）
- [x] T052 [P] ドキュメント更新（`specs/001-chem-model-webapp/quickstart.md`）
- [x] T053 [P] 開発補助: Justfile で Web/API 同時起動を追加（`Justfile`）
- [x] T054 [US1] .in パースの手動フォールバック追加（`apps/api/services/parse.py`）
- [x] T055 [P] フロント: 機能単位の分割と Shadcn/ui への移行（`apps/web/src/features/*`, `apps/web/src/components/ui/*`）
