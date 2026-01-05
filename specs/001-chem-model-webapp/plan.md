# Implementation Plan: 計算化学向け構造編集Webアプリ

**Branch**: `[001-chem-model-webapp]` | **Date**: 2026-01-05 | **Spec**: `specs/001-chem-model-webapp/spec.md`
**Input**: Feature specification from `specs/001-chem-model-webapp/spec.md`

## Summary

Quantum ESPRESSO の .in を中心に、構造の可視化・編集・比較・共有を行うWebアプリを構築する。フロントは TanStack Start + shadcn/ui + Mol*、バックエンドは FastAPI + ASE/pymatgen に分離する。P1（読み込み/編集/表示）を最優先で実装し、P2（並列比較・部分移植）とP3（スーパーセル・共有）を段階的に追加する。

## Technical Context

**Language/Version**: TypeScript (React), Python 3.11+  
**Primary Dependencies**: TanStack Start, shadcn/ui, Mol*, FastAPI, ASE, pymatgen  
**Storage**: 基本はファイル/メモリ（将来的に共有用ストレージ検討）  
**Testing**: Frontend: Vitest + Playwright(予定), Backend: pytest  
**Target Platform**: モダンブラウザ（Chrome/Firefox/Safari最新系） + ローカル/サーバ運用  
**Project Type**: Web (frontend + backend)  
**Performance Goals**: ~1k原子で 60fps 近い操作感、編集反映 <300ms  
**Constraints**: オフラインでも最低限の閲覧/編集が可能な設計、共有はリンク＋サーバ保存を優先し当面は単一HTML  
**Scale/Scope**: 2画面（Editor / Supercell）中心の中規模、構造数は任意

## Constitution Check

- `.specify/memory/constitution.md` が未整備のためゲートは未定義。必要なら先に合意済み原則を記載する。

## Project Structure

### Documentation (this feature)

```text
specs/001-chem-model-webapp/
├── spec.md
├── plan.md
├── research.md        # (必要に応じて追加)
├── data-model.md      # (必要に応じて追加)
├── quickstart.md      # (必要に応じて追加)
├── contracts/         # (API契約を追加)
└── tasks.md           # (Taskステージで作成)
```

### Source Code (repository root)

```text
apps/
├── web/               # TanStack Start + shadcn/ui + Mol*
└── api/               # FastAPI + ASE/pymatgen

packages/
└── shared/            # 型定義/ユーティリティ（必要なら）
```

**Structure Decision**: フロント/バック分離を前提に monorepo 構成。初期は apps/web と apps/api を用意し、共通型を packages/shared に配置する計画。

## Phased Plan

### Phase 0: 仕様確定（Spec承認）
- Spec のレビューと修正
- 共有方式（リンク＋サーバ保存優先、当面は単一HTML）を確定
- 取り扱う QE .in の範囲は当面「原子種・座標」に限定する

### Phase 1: 基盤設計
- データモデル設計（Structure/Atom/Lattice/Selection/ViewState）
- FastAPI API契約（parse, update, export, supercell）
- Mol* 表示統合設計（入力フォーマット変換、表示/色/ラベル/背景）

### Phase 2: MVP実装 (P1)
- .in 読み込み/パース（ASE/pymatgen）
- 3D表示 + テーブル編集の同期
- 基本編集操作（座標/元素変更、コピー/貼り付け）
- エクスポート（QE .in への書き戻し）

### Phase 3: 比較/移植 (P2)
- 複数構造の並列表示
- 部分選択コピー/移植
- 整列/シフト/距離計算

### Phase 4: スーパーセル/共有 (P3)
- 格子パラメータ入力/推定
- シーケンス指定によるスーパーセル生成
- 共有機能（単一リンク or 単一ファイル）

## Risks & Mitigations

- **QE .in の多様性**: ASE/pymatgen の対応範囲を明確化し、非対応は明示的エラーにする
- **Mol* への変換**: PDB/CIF 中間形式を採用し、変換処理を共通化
- **共有方式**: URL長やファイルサイズ上限に注意。大規模構造はサーバ保存方式を検討

## Open Questions

- 共有方式の実装（単一HTML → サーバ保存リンクへの移行順）
- 取り扱う QE .in の詳細仕様（constraints, magnetization, tags など）は後続
