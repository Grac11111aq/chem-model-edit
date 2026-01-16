# Feature Specification: ZPE 計算サーバ分離（Control Plane / Compute Plane）

**Feature Branch**: `[refactor/zpe-remote-compute]`  
**Created**: 2026-01-10  
**Status**: In Progress（HTTP仲介版に更新）  
**Input**: User description: "現状のZPE系の実装はローカル計算前提で、開発PC/サービス提供側で重い化学計算が走ってしまう。計算サーバー側とサービス提供側を分離したい。"

## ユーザーシナリオ & テスト（必須）

### User Story 1 - サービス提供側が計算せずにジョブを投入できる (Priority: P1)

UI から ZPE 計算を開始しても、サービス提供側ホストでは計算が走らず、計算サーバー側にジョブが委譲される。

**Why this priority**: 低スペック環境でのE2Eテスト、クラウドホスティング時の計算負荷回避が最重要のため。  
**Independent Test**: サービス提供側に QE/pseudo が存在しない状態でもジョブ投入が成功する。

**Acceptance Scenarios**:
1. **Given** サービス提供側に `pw.x` が無い、**When** `/calc/zpe/jobs` を呼ぶ、**Then** `job_id` が返る
2. **Given** 計算サーバーが有効、**When** ジョブ投入、**Then** 計算は計算サーバー側で実行される

---

### User Story 2 - 計算サーバーが結果を共有ストアへ保存する (Priority: P1)

計算サーバーが結果を共有ストアに保存し、サービス提供側はファイルシステムに依存せず結果取得できる。

**Why this priority**: 物理的に分離された環境間で結果を受け渡すため。  
**Independent Test**: サービス提供側に計算ディレクトリが無くても `/calc/zpe/jobs/{id}/result` が返る。

**Acceptance Scenarios**:
1. **Given** 計算完了ジョブ、**When** `/calc/zpe/jobs/{id}/result` を呼ぶ、**Then** 結果が JSON で返る
2. **Given** 計算完了ジョブ、**When** `/calc/zpe/jobs/{id}/files?kind=summary|freqs` を呼ぶ、**Then** CSV/summary が返る

---

### User Story 3 - E2Eテスト用に軽量なモック計算を選べる (Priority: P2)

計算が重い環境でもE2Eテストを回せるよう、モック計算モードを選択できる。

**Why this priority**: CI/E2E の実行環境に依存しない検証を可能にするため。  
**Independent Test**: モックモードで短時間に結果が返る。

**Acceptance Scenarios**:
1. **Given** モックモード、**When** `/calc/zpe/jobs` を呼ぶ、**Then** 数秒以内に完了ステータスが返る
2. **Given** モックモード、**When** 結果取得、**Then** 決定論的なダミー値が返る

---

### Edge Cases

- 計算サーバーがダウン/接続不可: MVP はジョブが queued のまま停滞する。結果取得は 409（未完了）を返し、詳細はログで追う。再試行/監視は後で追加する。
- 共有ストア（Redis/S3 など）障害: MVP は 503 を返して失敗を明確化し、例外はログに記録する。
- ワーカーの lease 期限切れ: 期限切れは再投入され、必要に応じて retry 上限で DLQ に送られる。
- 結果ファイルサイズが大きく共有ストア制限に抵触: MVP は Redis で扱える小容量前提。超過は失敗扱いとし、将来は S3 などへ移行する。
- ジョブの再実行（同一入力のキャッシュ）: MVP は明示的再投入のみ。キャッシュ戦略・再実行ポリシーは後で整理する。
- サービス提供側/計算サーバーで設定差分がある: MVP は結果に主要設定を埋め込み差分検出の材料にする。自動検知は後で追加する。

## 要件（必須）

### Functional Requirements

- **FR-001**: システムは `control-plane`（API）と `compute-plane`（計算）を論理的に分離する
- **FR-002**: `control-plane` は QE 実行系や pseudo を必要としない
- **FR-003**: 結果は共有ストアへ保存し、API はローカルファイルに依存しない
- **FR-004**: 計算バックエンドを `local` / `remote-http` / `remote-queue` / `mock` で切り替えられる
- **FR-005**: `remote-http` は control-plane の HTTPS API を介してワーカーがポーリングする
- **FR-006**: API は既存の `/calc/zpe/*` を維持し、UI 側の変更を最小化する
- **FR-007**: 通信はトークン認証または限定ネットワークで保護できる
- **FR-008**: 失敗時は原因（計算側/ストア側/設定不備）を明確に返す
- **FR-009**: 結果ストアは **control-plane のみが接続**し、ワーカーへ資格情報を配布しない
- **FR-010**: 計算サーバ登録は「短期登録トークン」を前提にし、登録後は **worker_token** で認証する  
  - 登録トークンは管理者が発行し、短期（例: 1時間）で失効する  
  - 配布は管理者が安全なチャネルで行い、期限切れトークンは無効扱い
- **FR-011**: ワーカーは lease を取得し、結果/失敗を HTTP で返却する
- **FR-012**: lease 期限切れ/失敗は retry 上限と DLQ を持つ  
  - 例: retry_max=3、指数バックオフ、上限超過は DLQ に送る（詳細は実装で定義）

### Key Entities

- **ZPEComputeBackend**: 計算委譲の切り替え責務を持つアダプタ
- **ResultStore**: 結果 JSON/CSV/summary を保存・取得するストア
- **ZPEJobRef**: job_id と結果参照情報（ストアキー等）
- **WorkerToken**: ワーカー認証に使うトークン
- **Lease**: job_id の一時的な実行権

## 成功条件（必須）

### Measurable Outcomes

- **SC-001**: サービス提供側は QE 無しで起動・ジョブ投入が可能
- **SC-002**: 計算結果取得は共有ストア経由で完結し、ローカルファイルに依存しない
- **SC-003**: モックモードのE2Eが 10 秒以内に完了する
- **SC-004**: 計算サーバー停止時に 5xx ではなく明確なエラー応答が返る
- **SC-005**: Upstash の URL/認証情報をユーザーへ配布しない

## E2E 検証手順（MVP）

詳細手順は `docs/zpe-worker-setup.ja.md` に集約し、Spec では要点のみを記載する。

1. **control-plane** で Redis と API を起動（`ZPE_COMPUTE_MODE=remote-http` / `ZPE_RESULT_STORE=redis`）
2. **admin token** を使って enroll token を発行し、計算サーバーを登録（worker_token を取得）
3. **compute-plane** で HTTP ワーカーを起動（`scripts/run-zpe-http-worker.sh`）
4. **H2 サンプル**（`samples/qe-in/h2_zpe.in`）を送信（`scripts/zpe-submit-h2.py`）
5. `/calc/zpe/jobs/{job_id}` をポーリングし、`/result` と `/files` が取得できることを確認

モック検証は `ZPE_COMPUTE_MODE=mock`（control-plane のみ）で行う。

## 決定事項

- 既存の RQ + Redis を維持しつつ、**HTTP仲介（remote-http）を主経路**とする
- 結果保存は初期は Redis へ格納（小さな CSV/summary 想定）し、将来 S3 などへ差し替え可能な設計にする
- API の `/calc/zpe/*` は継続し、UI 側を変更しない
- Redis は **control-plane のみが接続**し、ワーカーは HTTPS 経由でのみ接続する
- MVP では管理者が計算サーバ登録用の短期トークンを発行し、登録後は worker_token を利用する
