# ZPE リモートワーカー（HTTP仲介）Spec（ドラフト）

## 背景
- 現状は Redis(RQ) 直結で worker がジョブを取得・結果保存する構成。
- Upstash のURLは認証情報を含むため、ユーザーPCへ配布できない。
- ユーザーPCはNAT越しで動作し、Cloud Run から直接到達できない前提。

## 目的
- Upstashを開発者側（control-plane）にのみ保持し、**ワーカーはHTTPSでAPI経由**にジョブ取得/結果返却を行う。
- 現状のZPE API（フロントからのジョブ登録・結果取得）との整合性を維持する。

## 制約
- UpstashのURL/認証情報は**ユーザーに渡さない**。
- ユーザーPCは**外向きのHTTPS通信のみ**で動作する想定。
- 既存の `remote-queue`（Redis直結）モードは当面残す。
- H2など最小計算で動作検証できること。

## 想定構成（最小）
- **control-plane**: FastAPI（Cloud Run）
  - Upstash接続情報を保持
  - ジョブ登録APIを提供
  - ワーカー向けHTTPポーリングAPIを提供
- **compute-plane**: ユーザーPC（QE + worker）
  - API経由でジョブ取得
  - QE計算を実行
  - 結果をAPIへ返却
- **Redis**: Upstash（control-planeのみ接続）

## 仕様（案）

### 共通（HTTP）
- **認証**: `Authorization: Bearer <worker_token>`（worker向けAPI）
- **標準ステータス**:
  - `200 OK`: 正常
  - `204 No Content`: ジョブなし（lease）
  - `400 Bad Request`: 不正payload
  - `401 Unauthorized`: トークン不正/期限切れ
  - `403 Forbidden`: 権限不足/失効トークン
  - `404 Not Found`: job/workerが存在しない
  - `409 Conflict`: 既にリース済み/完了済み
  - `500 Internal Server Error`: サーバー側エラー
- **エラーレスポンス**:
```json
{
  "error": {
    "code": "INVALID_PAYLOAD",
    "message": "Missing required field: summary_text",
    "details": {}
  }
}
```

### Redisキー設計（案）
- `zpe:queue`: LIST（LPUSH/RPOP）
- `zpe:payload:{job_id}`: STRING（JSON）
- `zpe:status:{job_id}`: HASH（status/detail/updated_at）
- `zpe:lease:{job_id}`: HASH（worker_id/lease_id/expires_at）
- `zpe:lease:index`: ZSET（score=expires_at, member=job_id）
- `zpe:retry_count:{job_id}`: STRING（整数）
- `zpe:delay`: ZSET（score=run_at, member=job_id）
- `zpe:dlq`: LIST（失敗ジョブ）
- `zpe:result:{job_id}` / `zpe:summary:{job_id}` / `zpe:freqs:{job_id}`: STRING
- `zpe:token:{token_hash}`: HASH（worker_id/expires_at/revoked_at/label）
- `zpe:revoked_tokens`: SET（token_hash）
- TTL:
  - `zpe:lease:{job_id}`: **TTLは使わず** `zpe:lease:index` に期限を集約
  - `zpe:payload:{job_id}`: `result_ttl_seconds` と同等（完了後に削除も可）
  - `zpe:result/summary/freqs/status`: `result_ttl_seconds`

### 1. ジョブ投入（既存）
- `POST /calc/zpe/jobs`
- control-plane が `zpe:queue` に job_id を積む
- `zpe:payload:{job_id}` に payload（JSON）を保存
- `zpe:status:{job_id}` を `queued` に更新

### 2. ワーカー登録（拡張）
- 既存: `POST /calc/zpe/compute/enroll-tokens`
- 既存: `POST /calc/zpe/compute/servers/register`
- 追加案:
  - register時に **worker_token** を発行し、以後の認証に使用
  - **トークン形式**: まずは **opaque token**（ランダム文字列）を採用
  - **保存**: `zpe:token:{token_hash}` に `worker_id` を紐付けて保存
  - **有効期限**: 7日（`expires_at` を保存し、期限切れは 401）
  - **配布**: `/register` レスポンス本文に `worker_token` を含める
  - **バインド**: `worker_id` とトークンを必ず紐付け、検証時に一致チェック
  - **更新/ローテーション**: 期限切れ前に **再登録**（新しい enroll token で再register）
  - **失効**: `zpe:revoked_tokens` へ `token_hash` を登録（全リクエストでチェック）

### 3. ジョブ取得（新規）
- `POST /calc/zpe/compute/jobs/lease`
- 認証: `Bearer worker_token`
- 返却:
  - ジョブがある場合: `{ job_id, payload, lease_id, lease_ttl_seconds }`
  - ない場合: 204 No Content
- **lease取得は原子的に行う**:
  - Luaスクリプト等で `RPOP zpe:queue` → `HSET zpe:lease:{job_id}` を1回で実行
  - 成功時に `zpe:lease:index` に `expires_at` を登録（**TTLは使わない**）
- **lease_id**:
  - job取得時にランダム生成し `zpe:lease:{job_id}` に保存
  - **結果/失敗報告時の正当性検証**と **冪等性**に使用
- **ポーリング間隔**:
  - 推奨 5〜30秒（空の場合は指数バックオフ + jitter）

### 4. 結果返却（新規）
- `POST /calc/zpe/compute/jobs/{job_id}/result`
- 認証: `Bearer worker_token`
- payload:
```json
{
  "result": { "...": "..." },
  "summary_text": "string",
  "freqs_csv": "string",
  "meta": {
    "worker_hostname": "string",
    "computation_time_seconds": 12.3,
    "qe_version": "string",
    "timestamp": "ISO8601"
  }
}
```
- control-planeがUpstashへ保存:
  - **lease検証**: `zpe:lease:{job_id}` の worker_id/lease_id と一致すること
  - **冪等性**: 既に `status=finished` の場合は内容一致なら `200 OK`
  - **原子的書き込み**: result/summary/freqs/status をトランザクションまたはLuaで一括更新
  - **lease削除**: 成功時に `zpe:lease:{job_id}` と `zpe:lease:index` を削除

### 5. 失敗通知（新規）
- `POST /calc/zpe/compute/jobs/{job_id}/failed`
- payload:
```json
{
  "error_code": "QE_RUNTIME_ERROR",
  "error_message": "pw.x exited with code 1",
  "traceback": "optional"
}
```
- **lease検証**: `zpe:lease:{job_id}` の worker_id/lease_id と一致すること
- **失敗処理**:
  - retry回数を `zpe:retry_count:{job_id}` に記録
  - 上限未満なら `zpe:delay` に再投入（バックオフ）
  - 上限超過は `zpe:dlq` に移動し `status=failed` で確定
- **lease削除**: 失敗処理後に削除

### 6. 期限切れ処理（最小）
- **検知**: control-plane側で **定期リーパー**（例: 30秒周期）を実行
  - `zpe:lease:index` の期限切れを取得し、該当jobを再投入
  - `zpe:lease:{job_id}` が存在しない場合は **indexから除去してスキップ**
- **再投入**:
  - `zpe:retry_count:{job_id}` をインクリメント
  - **バックオフ**: `delay = min(300, 10 * 2^retry)` 秒
  - `zpe:delay` に `run_at` を保存し、到達時に `zpe:queue` へ移動
- **DLQ**:
  - retry上限（例: 3〜5回）超過は `zpe:dlq` へ移動

## セキュリティ方針
- Upstash URLはcontrol-planeのみ
- workerは **短期トークン** or **専用トークン**を使用
- トークンの失効（revocation）に備えてRedisで管理
- **失効API**: `DELETE /calc/zpe/compute/servers/{worker_id}` で該当トークンを失効
- **検証**: 全APIで `zpe:revoked_tokens` を確認

## 成功条件
- Upstash URLを**ユーザーに配布しない**
- ユーザーPCのworkerが **HTTPSでジョブ取得→QE計算→結果返却** できる
- フロント/APIから `GET /calc/zpe/jobs/{id}/result` で結果取得できる
- H2テストケースで ZPE が **0.26〜0.28 eV** 程度に収まる
- **Phase2**: ログ/メトリクス（queue depth, completion rate, active lease）を可視化

## 検証方法
- control-plane: FastAPI dev（この環境）
- compute-plane: 研究室PCで worker 起動
- Upstash: 本番相当の共有Redis
- H2計算でE2E確認
- **Phase2**: 複数worker同時実行、lease期限切れシナリオ、簡易負荷試験

## スコープ外（当面）
- AiiDA / Slurm 連携
- ジョブの可視化UI/詳細監視
- 大規模HPC最適化

## 次フェーズ（実装計画）
1. API仕様の確定（lease/timeout/retry）
2. control-planeのHTTP仲介API実装
3. workerのHTTPポーリング実装
4. docs/サンプル整備（H2入力）
5. E2E検証
