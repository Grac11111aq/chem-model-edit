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

### 1. ジョブ投入（既存）
- `POST /calc/zpe/jobs`
- control-plane が `zpe:queue` に job_id を積む
- `zpe:payload:{job_id}` に payload を保存
- `zpe:status:{job_id}` を `queued` に更新

### 2. ワーカー登録（拡張）
- 既存: `POST /calc/zpe/compute/enroll-tokens`
- 既存: `POST /calc/zpe/compute/servers/register`
- 追加案:
  - register時に **worker_token** を発行し、以後の認証に使用
  - `Authorization: Bearer <worker_token>` 必須

### 3. ジョブ取得（新規）
- `POST /calc/zpe/compute/jobs/lease`
- 認証: `Bearer worker_token`
- 返却:
  - ジョブがある場合: `{ job_id, payload, lease_id, lease_ttl_seconds }`
  - ない場合: 204 No Content
- Redis側:
  - `zpe:lease:{job_id}` に worker_id + TTL を保存
  - leaseがある場合は他workerに渡さない

### 4. 結果返却（新規）
- `POST /calc/zpe/compute/jobs/{job_id}/result`
- 認証: `Bearer worker_token`
- payload: `{ result, summary_text, freqs_csv, meta? }`
- control-planeがUpstashへ保存:
  - `zpe:result:{job_id}`
  - `zpe:summary:{job_id}`
  - `zpe:freqs:{job_id}`
  - `zpe:status:{job_id} = finished`

### 5. 失敗通知（新規）
- `POST /calc/zpe/compute/jobs/{job_id}/failed`
- 失敗理由を保存し `status=failed`

### 6. 期限切れ処理（最小）
- lease TTLが切れた場合は再度 `queue` に戻す（簡易）
- 最初は **シンプルな再投入** を採用

## セキュリティ方針
- Upstash URLはcontrol-planeのみ
- workerは **短期トークン** or **専用トークン**を使用
- トークンの失効（revocation）に備えてRedisで管理

## 成功条件
- Upstash URLを**ユーザーに配布しない**
- ユーザーPCのworkerが **HTTPSでジョブ取得→QE計算→結果返却** できる
- フロント/APIから `GET /calc/zpe/jobs/{id}/result` で結果取得できる
- H2テストケースで ZPE が **0.26〜0.28 eV** 程度に収まる

## 検証方法
- control-plane: FastAPI dev（この環境）
- compute-plane: 研究室PCで worker 起動
- Upstash: 本番相当の共有Redis
- H2計算でE2E確認

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
