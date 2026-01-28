# ZPE ワーカーセットアップ（compute-plane）

このガイドは、ZPE 計算ワーカーを別マシンで動かすための手順です。
control-plane（FastAPI）はジョブ投入のみを行い、QE 計算は compute-plane に委譲します。

## 構成
- **control-plane**: FastAPI API（Cloud Run など）
- **compute-plane**: HTTP worker（別マシン）
- **共有ストア**: Redis（control-plane のみが接続）

## 秘密情報の境界
- **Admin token** は control-plane のみで保持
- **compute-plane** は短期 enroll token を使って登録

## 前提（compute-plane）
- Python >= 3.13（`uv` 使用）
- Quantum ESPRESSO（`pw.x`）が利用可能
- Pseudopotential ディレクトリ
- control-plane へ HTTPS で到達可能

## リポジトリと Python 環境
```bash
# SSH (事前にSSH鍵の設定が必要)
git clone git@github.com:grace-bidos/chem-model-edit.git
# or HTTPS
# git clone https://github.com/grace-bidos/chem-model-edit.git
cd chem-model-edit/apps/api
uv sync
```

## 環境変数ファイル
control-plane と compute-plane を分離します。

- control-plane: `apps/api/.env.control.example`
- compute-plane: `apps/api/.env.compute.example`

### Control-plane（FastAPI）
remote-http 用の最小例:
```bash
ZPE_REDIS_URL=redis://localhost:6379/0
ZPE_QUEUE_NAME=zpe
ZPE_COMPUTE_MODE=remote-http
ZPE_RESULT_STORE=redis
ZPE_ADMIN_TOKEN=change-me
```

### Compute-plane（worker）
QE 実行に必要な最小例:
```bash
ZPE_CONTROL_API_URL=http://localhost:8000
ZPE_WORKER_TOKEN=change-me
ZPE_PSEUDO_DIR=/path/to/pseudo
ZPE_PW_PATH=/path/to/pw.x
```

任意設定:
- `ZPE_USE_MPI=true|false`
- `ZPE_MPI_CMD=mpirun`
- `ZPE_NP_CORE=12`
- `ZPE_WORK_DIR=~/zpe_jobs`
- `ZPE_ALLOW_INPUT_PSEUDO_DIR=false`
- `ZPE_ENVIRON_PATH=/path/to/environ.in`
- `ZPE_WORKER_MODE=mock`（QE を回避してダミー結果）

## Enroll token フロー（任意だが推奨）

### ユーザー自己発行（UI経由）
1) Web UI でサインインし、Enroll token を発行（有効期限は約 1 時間）。
2) compute-plane 登録時に queue 名を指定して登録:
```bash
curl -X POST http://localhost:8000/calc/zpe/compute/servers/register \
  -H "Content-Type: application/json" \
  -d '{"token": "<ENROLL_TOKEN>", "name": "worker-1", "queue_name": "zpe", "meta": {"host": "compute-01"}}'
```
登録された worker は UI のキューターゲット一覧に表示されます。

### 1) control-plane でトークン発行
```bash
curl -X POST http://localhost:8000/calc/zpe/compute/enroll-tokens \
  -H "Content-Type: application/json" \
  -H "Authorization: Bearer $ZPE_ADMIN_TOKEN" \
  -d '{"ttl_seconds": 3600, "label": "worker-1"}'
```

### 2) compute-plane 登録
```bash
curl -X POST http://localhost:8000/calc/zpe/compute/servers/register \
  -H "Content-Type: application/json" \
  -d '{"token": "<ENROLL_TOKEN>", "name": "worker-1", "queue_name": "zpe", "meta": {"host": "compute-01"}}'
```

## Worker 起動
補助スクリプトを用意しています:
```bash
./scripts/run-zpe-worker.sh
```

HTTP ワーカーの場合:
```bash
./scripts/run-zpe-http-worker.sh
```

特定の env ファイルを指定する場合:
```bash
ZPE_ENV_FILE=apps/api/.env.compute ./scripts/run-zpe-worker.sh
```

`just` を使う場合:
```bash
just zpe-worker
```

HTTP ワーカーの場合:
```bash
just zpe-http-worker
```

起動時に `uv sync` も実行する場合:
```bash
ZPE_WORKER_SYNC=1 ./scripts/run-zpe-worker.sh
```

## スモークテスト（E2E）
1) API 起動（control-plane）
```bash
cd apps/api
uv run uvicorn main:app --reload --port 8000
```

2) ジョブ投入（別ターミナル・repo root で実行）
```bash
python - <<'PY'
import json
import urllib.request
from pathlib import Path

content = Path("samples/qe-in/Al001_m4_relax_fromscratch.in").read_text()
payload = {
    "content": content,
    "mobile_indices": [0],
    "calc_mode": "new",
    "use_environ": False,
}
req = urllib.request.Request(
    "http://localhost:8000/calc/zpe/jobs",
    data=json.dumps(payload).encode(),
    headers={
        "Content-Type": "application/json",
        "Authorization": "Bearer <USER_SESSION_TOKEN>",
    },
)
print(urllib.request.urlopen(req).read().decode())
PY
```

3) ステータス確認と結果取得
```bash
curl http://localhost:8000/calc/zpe/jobs/<JOB_ID>
curl http://localhost:8000/calc/zpe/jobs/<JOB_ID>/result
```

## Mock mode（control-plane のみ）
E2E/CI 向けに worker を不要にする場合:
```bash
ZPE_COMPUTE_MODE=mock
```
決定論的なダミー結果が返るため、compute-plane は不要です。

## Worker 側の mock 実行（遠隔経路の検証用）
**遠隔キュー経路**を通しつつ QE を避けたい場合は compute-plane に以下を設定します:
```bash
ZPE_WORKER_MODE=mock
```
API 側は `ZPE_COMPUTE_MODE=remote-queue` のままにします。
