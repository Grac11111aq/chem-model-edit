# Implementation Plan: 計算化学向け構造編集Webアプリ

**Branch**: `[001-chem-model-webapp]` | **Date**: 2026-01-05 | **Spec**: `specs/001-chem-model-webapp/spec.md`
**Input**: Feature specification from `specs/001-chem-model-webapp/spec.md`

## Summary

Quantum ESPRESSO の .in を中心に、2構造(A/B)の可視化・編集・表面転写・複合スーパーセル構築・距離QC・格子変換・共有/重ね表示を行うWebアプリを構築する。フロントは TanStack Start + shadcn/ui + Mol*、バックエンドは FastAPI + ASE/pymatgen に分離する。P1（読み込み/編集/表示）を最優先で実装し、P2（表面転写/距離QC）とP3（タイル合成/格子変換/共有）を段階的に追加する。

## Technical Context

**Language/Version**: TypeScript (React), Python 3.11+  
**Primary Dependencies**: TanStack Start, shadcn/ui, Mol*, FastAPI, ASE, pymatgen  
**Storage**: 基本はファイル/メモリ  
**Testing**: Frontend: Vitest + Playwright(予定), Backend: pytest  
**Target Platform**: モダンブラウザ（Chrome/Firefox/Safari最新系） + ローカル/サーバ運用  
**Project Type**: Web (frontend + backend)  
**Performance Goals**: ~1k原子で 60fps 近い操作感、編集反映 <300ms  
**Constraints**: 対象は QE .in、出力は入力形式/単位を保持がデフォルト、celldmは保持（再計算/削除は将来対応）、共有は当面単一HTML  
**Scale/Scope**: 2画面（Editor / Supercell）中心の中規模、複数構造IDの扱いを前提（Editor v2 対応）

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
- 出力オプションは設定画面で管理し、デフォルトは入力形式/単位保持を確定
- celldm は保持のみ（再計算/削除は将来対応）を明記
- 取り扱う QE .in の範囲は当面「原子種・座標・格子」に限定する
- 共有は単一HTML（重ね表示/透明度/可視を再現）を確定

### Phase 1: 基盤設計
- データモデル設計（Model(A/B)/Structure/Atom/Lattice/Selection/SurfaceTransfer/DistanceReport/TilingPattern/ExportSettings）
- StructureRegistry と StructureId ベースの契約（Editor v2/Supercell v2）
- FastAPI API契約（parse, export, transfer, distance, supercell, lattice-convert）
- Mol* 表示統合設計（入力フォーマット変換、表示/色/ラベル/背景）

### Phase 2: MVP実装 (P1)
- .in 読み込み/パース（ASE/pymatgen）
- 3D表示 + テーブル編集の同期
- 基本編集操作（座標/元素変更、コピー/貼り付け）
- エクスポート（QE .in への書き戻し、入力形式/単位保持）
- 出力設定画面（デフォルト保持のみ、将来拡張を見据えたUI枠）

### Phase 2b: Editor UI刷新（/editor-v2 試行）
- 新UIモックを `/editor-v2` に実装し、既存 `/editor` は維持
- 左ツールナビ + ファイルマネージャ + ワークスペース構成へ移行
- 特別モードの命名を整理（内部ID: transfer/supercell/vibration、表示名: Transfer / Supercell Builder / Vibrations (Preview)）
- 機能接続は段階的に進行（最初は見た目・遷移、次にMol*・編集系）
- ZPE/振動系は当面プレースホルダで接続なし

### Phase 3: 表面転写/距離QC (P2)
- 2構造(A/B)の並列表示
- 選択UI（表面領域の選択を含む）
- インデックス対応の表面転写（Cartesian座標上書き）
- 対応原子距離の一覧（PBC時は最小像）
- A/B不整合時の事前警告（原子数・格子非互換）
- 重ね表示モードの切り替えと構造ごとの表示/透明度調整

### Phase 3b: Δ移植（小スラブ out → 大スラブ in）(P2)
- API追加: `POST /transplant/delta`
  - 入力: `small_in`, `small_out`, `large_in`（全て文字列）
  - 出力: `content`（移植済み .in）
  - 仕様: 小スラブ出力の初期/最終座標で変位を算出し、初期座標が一致する大スラブ原子（元素一致）へ適用
  - エラー: 可動フラグ欠落/初期座標の一致失敗（未発見・重複）/ATOMIC_POSITIONS欠落
  - 制約: 非直交セル未対応（wrap は Lx/Ly のみ）
- UI追加: 専用ページ（例: `/transplant`）
  - 3入力（Small .in / Small .out / Large .in）
  - 実行ボタン、結果プレビュー、コピー/ダウンロード
  - エラー表示（APIエラーの内容を明示）
- テスト: APIユニット + APIルートテスト（正常/異常系）
- Editor v2: 移植ツールで .out を入力し Δ移植を実行
  - Source/Target の .in（qeInput）と .out を組み合わせて API を呼び出す
  - 結果プレビューと .in ダウンロードを接続

### Phase 4: 複合スーパーセル/格子変換 (P3)
- タイルパターン入力（StructureId の 2D グリッド or GUI）
- タイル合成（タイル内座標保持＋並進のみ、格子は baseStructure の lattice を利用）
- 重複/衝突チェック（許容誤差つき、任意）
- 格子パラメータ変換（格子定数 ⇄ セルベクトル）
- 共有HTML出力（重ね表示/可視/透明度の再現）

## Risks & Mitigations

- **QE .in の多様性**: ASE/pymatgen の対応範囲を明確化し、非対応は明示的エラーにする
- **Mol* への変換**: PDB/CIF 中間形式を採用し、変換処理を共通化
- **インデックス対応の前提**: A/Bの整列・並び順が崩れると誤転写の恐れがあるため、事前警告とガイドを用意
- **celldm 整合性**: 出力保持の前提で、セル編集時の不整合に注意喚起を出す
- **共有HTMLの容量**: 大規模構造ではサイズが肥大化するため、将来はサーバ保存方式を検討

## Open Questions

- タイルパターン入力のUI設計（2D配列/GUIの詳細）
- 表面選択の補助（Z閾値など）を導入するか
- celldm の再計算/削除をどの段階で実装するか
- 取り扱う QE .in の詳細仕様（constraints, magnetization, tags など）は後続
- 共有のサーバ保存リンク方式は将来検討
