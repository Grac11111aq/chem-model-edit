# Feature Specification: ZPEGUI バックエンド移植（ZPE/Vibrations）

**Feature Branch**: `[dev]`  
**Created**: 2026-01-07  
**Status**: Draft  
**Input**: User description: "ref-legacy/ZPEGUI9.py の機能をバックエンドへ移植（FastAPI + RQ + ASE/pymatgen）。一旦バックエンドのみ実装。"

## ユーザーシナリオ & テスト（必須）

### User Story 1 - QE入力を解析し可動/固定原子情報を返す (Priority: P1)

ユーザーが QE の .in をアップロードすると、構造（原子/格子）と固定原子の情報が得られる。

**Why this priority**: 以降の可動原子選択と計算に必須の前処理であるため。  
**Independent Test**: QE .in を投げ、原子数・座標・固定フラグが期待通りに返る。

**Acceptance Scenarios**:
1. **Given** 有効な QE .in、**When** `/calc/zpe/parse` を呼ぶ、**Then** 原子一覧・格子・固定原子 indices が返る
2. **Given** ATOMIC_POSITIONS が欠落した .in、**When** `/calc/zpe/parse` を呼ぶ、**Then** 400 でエラーが返る

---

### User Story 2 - ZPE計算ジョブを投入し、非同期で計算できる (Priority: P1)

ユーザーが可動原子 indices を指定して計算を開始すると、API がジョブをキューに積み、ワーカーが計算を実行できる。

**Why this priority**: ZPE/Vibrations 実行が主要機能であり、UIは非同期前提のため。  
**Independent Test**: ジョブ投入時に即時 job_id が返り、RQ にエンキューされる。

**Acceptance Scenarios**:
1. **Given** 構造と可動原子 indices、**When** `/calc/zpe/jobs` を呼ぶ、**Then** job_id が返る
2. **Given** 可動原子が空、**When** `/calc/zpe/jobs` を呼ぶ、**Then** 400 でエラーが返る
3. **Given** qe/pseudo の必須ファイルが不足、**When** ジョブ実行、**Then** job は failed になり理由が取得できる

---

### User Story 3 - ジョブの状態/結果を取得し、CSVも入手できる (Priority: P2)

ユーザーがジョブの状態をポーリングし、完了後に周波数・ZPE・S_vib と CSV を取得できる。

**Why this priority**: UIの表示（グラフ/CSV DL）に必要なため。  
**Independent Test**: 完了ジョブの result が JSON で返り、CSV 取得エンドポイントが機能する。

**Acceptance Scenarios**:
1. **Given** 実行中ジョブ、**When** `/calc/zpe/jobs/{id}` を呼ぶ、**Then** status=started/queued が返る
2. **Given** 完了ジョブ、**When** `/calc/zpe/jobs/{id}/result` を呼ぶ、**Then** ZPE/周波数/メタ情報が返る
3. **Given** 完了ジョブ、**When** `/calc/zpe/jobs/{id}/files?kind=summary|freqs` を呼ぶ、**Then** summary.txt / freqs.csv がDLできる

---

### Edge Cases

- `pw.x` が見つからない/実行できない
- pseudopotentials が不足・不一致
- ATOMIC_POSITIONS の末尾フラグが不正（固定判定不能）
- ジョブの作業ディレクトリ破損/途中停止
- Redis が停止・接続不可
- ジョブ再実行（new/continue）の扱い

## 要件（必須）

### Functional Requirements

- **FR-001**: システムは QE .in を解析し、原子・格子・固定原子 indices を返す
- **FR-002**: システムは可動原子 indices を受け取り、ZPE計算ジョブを RQ に投入できる
- **FR-003**: ワーカーは ASE Vibrations + QE を用いて周波数/ZPE/S_vib を算出する
- **FR-004**: ジョブは作業ディレクトリ単位で隔離し、キャッシュ健全性チェックを行う
- **FR-005**: ジョブ状態/結果取得APIを提供し、エラー時は理由を返す
- **FR-006**: 結果は JSON と CSV/summary で取得できる
- **FR-007**: environ.in 使用の有無を指定でき、必要な場所へ配布できる
- **FR-008**: QE入力から system/electrons/kpoints(automatic) を読み取り、計算に反映できる
- **FR-009**: pseudopotentials は `.env` 指定のディレクトリを優先し、許可時のみ入力の pseudo_dir を解決できる

### Key Entities

- **ZPEJob**: 解析対象構造・可動原子・設定・作業ディレクトリを持つジョブ
- **ZPEJobResult**: 周波数配列、ZPE、S_vib、実行時間、使用設定
- **ZPEParseResult**: 構造データ + 固定原子 indices

## 成功条件（必須）

### Measurable Outcomes

- **SC-001**: 1k 原子の .in を 2 秒以内に解析できる
- **SC-002**: ジョブ投入 API が 500ms 以内に job_id を返す
- **SC-003**: 完了ジョブが ZPE/周波数/メタ情報を JSON で取得できる
- **SC-004**: 失敗ジョブで原因が明確に返る（例: pw.x 不在）

## 決定事項

- pseudopotentials は `.env` 指定のディレクトリを利用（必要なら input の pseudo_dir を許可）
- `use_mpi` / `np_core` は `.env` で固定（ローカル利用前提）
- 保存先は Redis + job 作業ディレクトリ（SQLiteは使用しない）
- エンドポイントは `calc/zpe` プレフィックスで統一
