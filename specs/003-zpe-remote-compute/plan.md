# Implementation Plan: ZPE 計算サーバ分離

**Branch**: `[refactor/zpe-remote-compute]` | **Date**: 2026-01-10 | **Spec**: `specs/003-zpe-remote-compute/spec.md`
**Input**: Feature specification from `specs/003-zpe-remote-compute/spec.md`

## Summary

ZPE 計算を `control-plane`（API）と `compute-plane`（計算）に分離し、サービス提供側では計算が走らない構成にする。ジョブは **HTTP 仲介（remote-http）** でワーカーがポーリングし、結果は control-plane が共有ストア（初期は Redis）へ保存する。E2E/CI 向けにモック計算モードも提供する。

## Technical Context

- **Language/Version**: Python 3.13
- **Primary Dependencies**: FastAPI, RQ, redis-py, pydantic-settings
- **Result Store**: control-plane のみが接続する Redis（将来 S3 へ差し替え可能な設計）
- **Testing**: pytest（モックバックエンド + ストアのユニットテスト）
- **Target Platform**: API: 任意のクラウド / Compute: 計算サーバー
- **Constraints**: API 側は `pw.x` や pseudo を持たない前提

## Project Structure

### Documentation (this feature)

```plaintext
specs/003-zpe-remote-compute/
├── spec.md
├── plan.md
└── tasks.md
```

### Source Code (repository root)

```plaintext
apps/api/
├── services/zpe/
│   ├── backends.py      # 計算委譲アダプタ
│   ├── result_store.py  # 結果保存/取得
│   └── settings.py      # compute mode など追加
└── main.py              # /calc/zpe/* でバックエンド切替を利用
```

**Structure Decision**: 既存の `apps/api/services/zpe` にアダプタ層を追加し、API 自体は `/calc/zpe/*` を維持する。

## Work Phases

1. **Design/Setup**
   - `ZPE_COMPUTE_MODE` / `ZPE_RESULT_STORE` など設定追加
   - Public Redis + TLS + ACL 前提の構成を明示
   - 短期登録トークン方式（将来のユーザー単位登録へ拡張）を設計
   - 結果ストアとバックエンドのインターフェース定義
2. **Core Implementation**
   - Job enqueue API（`remote-http` backend）
   - ワーカー認証（worker_token）と登録/失効 API
   - Lease 取得 API と期限切れ処理
   - 結果返却/失敗返却 API
   - Retry/backoff/DLQ の最小実装
   - HTTP ワーカーポーリング実装
   - `remote-queue` / `mock` バックエンド（互換性維持）
   - 結果保存/取得の統一（Redis）
3. **Verification**
   - HTTP 仲介経由の結果取得テスト
   - モックモードのE2E向け検証

## Observability (MVP)

- control-plane は job_id とステータス遷移、RedisError をログに残す
- compute-plane は job_id を含む開始/完了/失敗ログを出す
- 監視ダッシュボード・アラート・SLO は運用が安定してから追加する
- **Structured Logging（Planned）**
  - JSON 形式の構造化ログを採用し、control-plane / compute-plane で相関IDを揃える
  - 最低限のフィールド: request_id, job_id, calc_id, stage, status, exit_code, duration_ms, qe_version, backend/result_store, user_id（任意）
  - 実装は Issue #47 で追跡する

## Rollback Plan

- `ZPE_COMPUTE_MODE` を `local` に戻す
- compute-plane ワーカーを停止する
- `/calc/zpe/*` の挙動を従来のローカル方式に戻す
- Redis 側のジョブ/結果はそのまま残す（必要に応じて手動削除）
- **in-flight ジョブ**: 新規受付停止 → 実行中ジョブは完了待ち（または強制停止してDLQへ退避）
- 段階的切替やデータクリーンアップの詳細手順は運用が固まってから追記する
- **運用ガードレール（Planned）**
  - 新規 remote 受付を止めるフラグ/トグル（既存ジョブは継続 or 明示キャンセル）
  - worker の dequeue を止める運用手順
  - 実装は Issue #47 で追跡する
