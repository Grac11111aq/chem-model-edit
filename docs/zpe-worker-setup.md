# ZPE Worker Setup (compute-plane)

This guide explains how to run the ZPE compute worker on a separate machine.
The control-plane (FastAPI) only enqueues jobs and never runs QE locally.

## Architecture
- **Control-plane**: FastAPI API on Cloud Run (or any host)
- **Compute-plane**: HTTP worker on a separate machine
- **Shared store**: Redis (control-plane only)

## Secret boundary
- **Admin token** stays in the control-plane only.
- **Compute-plane** uses a short-lived enroll token to register itself.

## Prerequisites (compute-plane)
- Python >= 3.13 with `uv`
- Quantum ESPRESSO (`pw.x`) installed
- Pseudopotential directory
- Reachable control-plane over HTTPS

## Repo + Python env
```bash
git clone git@github.com:grace-bidos/chem-model-edit.git
# or HTTPS
# git clone https://github.com/grace-bidos/chem-model-edit.git
cd chem-model-edit/apps/api
uv sync
```

## Environment files
Use separate env files for control-plane and compute-plane.

- Control-plane example: `apps/api/.env.control.example`
- Compute-plane example: `apps/api/.env.compute.example`

### Control-plane (FastAPI)
Minimal config for remote-http:
```bash
ZPE_REDIS_URL=redis://localhost:6379/0
ZPE_QUEUE_NAME=zpe
ZPE_COMPUTE_MODE=remote-http
ZPE_RESULT_STORE=redis
ZPE_ADMIN_TOKEN=change-me
```

### Compute-plane (worker)
Minimal config for QE execution:
```bash
ZPE_CONTROL_API_URL=http://localhost:8000
ZPE_WORKER_TOKEN=change-me
ZPE_PSEUDO_DIR=/path/to/pseudo
ZPE_PW_PATH=/path/to/pw.x
```

Optional:
- `ZPE_USE_MPI=true|false`
- `ZPE_MPI_CMD=mpirun`
- `ZPE_NP_CORE=12`
- `ZPE_WORK_DIR=~/zpe_jobs`
- `ZPE_ALLOW_INPUT_PSEUDO_DIR=false`
- `ZPE_ENVIRON_PATH=/path/to/environ.in`
- `ZPE_WORKER_MODE=mock` (skip QE; use dummy results)

## Enroll-token flow (optional but recommended)

### User self-service (UI-driven)
1) Sign in on the web UI and generate an enroll token (valid for ~1 hour).
2) Use the token on the compute server registration call with a queue name:
```bash
curl -X POST http://localhost:8000/calc/zpe/compute/servers/register \
  -H "Content-Type: application/json" \
  -d '{"token": "<ENROLL_TOKEN>", "name": "worker-1", "queue_name": "zpe", "meta": {"host": "compute-01"}}'
```
The server is registered under the signed-in user and becomes available as a queue target in the UI.

### 1) Create a token on the control-plane
```bash
curl -X POST http://localhost:8000/calc/zpe/compute/enroll-tokens \
  -H "Content-Type: application/json" \
  -H "Authorization: Bearer $ZPE_ADMIN_TOKEN" \
  -d '{"ttl_seconds": 3600, "label": "worker-1"}'
```

### 2) Register the compute server
```bash
curl -X POST http://localhost:8000/calc/zpe/compute/servers/register \
  -H "Content-Type: application/json" \
  -d '{"token": "<ENROLL_TOKEN>", "name": "worker-1", "queue_name": "zpe", "meta": {"host": "compute-01"}}'
```

## Start the worker
A helper script is provided:
```bash
./scripts/run-zpe-worker.sh
```

For HTTP worker:
```bash
./scripts/run-zpe-http-worker.sh
```

You can point it to a specific env file:
```bash
ZPE_ENV_FILE=apps/api/.env.compute ./scripts/run-zpe-worker.sh
```

If you use `just`:
```bash
just zpe-worker
```

For HTTP worker:
```bash
just zpe-http-worker
```

If you need to run `uv sync` at startup, set:
```bash
ZPE_WORKER_SYNC=1 ./scripts/run-zpe-worker.sh
```

## Smoke test (end-to-end)
1) Start the API (control-plane)
```bash
cd apps/api
uv run uvicorn main:app --reload --port 8000
```

2) Submit a job (in another terminal, repo root)
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

3) Poll the job and fetch the result
```bash
curl http://localhost:8000/calc/zpe/jobs/<JOB_ID>
curl http://localhost:8000/calc/zpe/jobs/<JOB_ID>/result
```

## Mock mode (control-plane only)
For fast E2E/CI without a worker:
```bash
ZPE_COMPUTE_MODE=mock
```
The API returns deterministic dummy results, and no worker is required.

## Mock execution on the worker (remote-path test)
If you want to test the **remote queue path** without running QE, set this on the compute-plane:
```bash
ZPE_WORKER_MODE=mock
```
The worker will dequeue jobs and write dummy results to Redis, while the API remains
in `ZPE_COMPUTE_MODE=remote-queue`.
