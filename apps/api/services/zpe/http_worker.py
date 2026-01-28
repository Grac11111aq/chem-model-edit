from __future__ import annotations

import json
import socket
import time
import traceback
from datetime import datetime, timezone
from typing import Any, Dict, Optional, Tuple
from urllib import request as urlrequest
from urllib import error as urlerror

from .settings import get_zpe_settings
from .worker import compute_zpe_artifacts


def _now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def _request_json(
    method: str,
    url: str,
    *,
    token: str,
    payload: Optional[Dict[str, Any]] = None,
    timeout: int = 60,
) -> Tuple[int, Optional[Dict[str, Any]]]:
    headers = {
        "Authorization": f"Bearer {token}",
        "Accept": "application/json",
    }
    data = None
    if payload is not None:
        headers["Content-Type"] = "application/json"
        data = json.dumps(payload).encode("utf-8")
    req = urlrequest.Request(url, data=data, headers=headers, method=method)
    try:
        with urlrequest.urlopen(req, timeout=timeout) as resp:
            status = resp.status
            body = resp.read().decode("utf-8")
    except urlerror.HTTPError as exc:
        status = exc.code
        body = exc.read().decode("utf-8") if exc.fp else ""
    except urlerror.URLError:
        return 0, None
    if not body:
        return status, None
    try:
        return status, json.loads(body)
    except json.JSONDecodeError:
        return status, None


def _lease_job(base_url: str, token: str, timeout: int) -> Optional[Dict[str, Any]]:
    status, data = _request_json(
        "POST",
        f"{base_url}/calc/zpe/compute/jobs/lease",
        token=token,
        payload=None,
        timeout=timeout,
    )
    if status == 204:
        return None
    if status != 200:
        raise RuntimeError(f"lease failed (status={status})")
    return data


def _submit_result(
    base_url: str,
    token: str,
    job_id: str,
    lease_id: str,
    result: Dict[str, Any],
    summary_text: str,
    freqs_csv: str,
    meta: Dict[str, Any],
    timeout: int,
) -> None:
    payload = {
        "lease_id": lease_id,
        "result": result,
        "summary_text": summary_text,
        "freqs_csv": freqs_csv,
        "meta": meta,
    }
    status, _ = _request_json(
        "POST",
        f"{base_url}/calc/zpe/compute/jobs/{job_id}/result",
        token=token,
        payload=payload,
        timeout=timeout,
    )
    if status not in {200, 409}:
        raise RuntimeError(f"result submit failed (status={status})")


def _submit_failed(
    base_url: str,
    token: str,
    job_id: str,
    lease_id: str,
    error_code: str,
    error_message: str,
    tb: str,
    timeout: int,
) -> None:
    payload = {
        "lease_id": lease_id,
        "error_code": error_code,
        "error_message": error_message,
        "traceback": tb,
    }
    _request_json(
        "POST",
        f"{base_url}/calc/zpe/compute/jobs/{job_id}/failed",
        token=token,
        payload=payload,
        timeout=timeout,
    )


def run_http_worker() -> None:
    settings = get_zpe_settings()
    if not settings.worker_token:
        raise RuntimeError("ZPE_WORKER_TOKEN is required for HTTP worker")

    base_url = settings.control_api_url.rstrip("/")
    token = settings.worker_token
    poll_interval = max(1, int(settings.worker_poll_interval_seconds))
    max_interval = max(poll_interval, int(settings.worker_poll_max_interval_seconds))
    timeout = max(1, int(settings.worker_request_timeout_seconds))

    backoff = poll_interval
    hostname = socket.gethostname()

    while True:
        try:
            lease = _lease_job(base_url, token, timeout)
        except Exception as exc:
            print(f"[zpe-http-worker] lease error: {exc}")
            time.sleep(backoff)
            backoff = min(max_interval, backoff * 2)
            continue

        if lease is None:
            time.sleep(backoff)
            backoff = min(max_interval, backoff * 2)
            continue

        backoff = poll_interval
        job_id = lease.get("job_id")
        payload = lease.get("payload")
        lease_id = lease.get("lease_id")
        if not job_id or not payload or not lease_id:
            print("[zpe-http-worker] invalid lease payload")
            time.sleep(poll_interval)
            continue

        start = time.monotonic()
        try:
            artifacts = compute_zpe_artifacts(payload, job_id=job_id)
            elapsed = time.monotonic() - start
            meta = {
                "worker_hostname": hostname,
                "computation_time_seconds": elapsed,
                "timestamp": _now_iso(),
            }
            _submit_result(
                base_url,
                token,
                job_id,
                lease_id,
                artifacts.result,
                artifacts.summary_text,
                artifacts.freqs_csv,
                meta,
                timeout,
            )
        except Exception as exc:
            tb = traceback.format_exc()
            print(f"[zpe-http-worker] job failed: {exc}")
            _submit_failed(
                base_url,
                token,
                job_id,
                lease_id,
                "COMPUTE_ERROR",
                str(exc),
                tb,
                timeout,
            )


if __name__ == "__main__":
    run_http_worker()
