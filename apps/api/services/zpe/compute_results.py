from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
import json
from typing import Any, Dict, Optional

from redis import Redis, WatchError

from .queue import get_redis_connection
from .settings import get_zpe_settings


_STATUS_PREFIX = "zpe:status:"
_RESULT_PREFIX = "zpe:result:"
_SUMMARY_PREFIX = "zpe:summary:"
_FREQS_PREFIX = "zpe:freqs:"
_LEASE_PREFIX = "zpe:lease:"
_LEASE_INDEX = "zpe:lease:index"
_RETRY_PREFIX = "zpe:retry_count:"
_DELAY_ZSET = "zpe:delay"
_DLQ = "zpe:dlq"


def _now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def _decode_map(data: Dict[bytes, bytes]) -> Dict[str, str]:
    return {k.decode("utf-8"): v.decode("utf-8") for k, v in data.items()}


def _get_lease(redis: Redis, job_id: str) -> Dict[str, str]:
    raw = redis.hgetall(f"{_LEASE_PREFIX}{job_id}")
    return _decode_map(raw) if raw else {}


@dataclass
class ResultSubmitOutcome:
    idempotent: bool


@dataclass
class FailureOutcome:
    requeued: bool
    retry_count: int


def submit_result(
    *,
    job_id: str,
    worker_id: str,
    lease_id: str,
    result: Dict[str, Any],
    summary_text: str,
    freqs_csv: str,
) -> ResultSubmitOutcome:
    settings = get_zpe_settings()
    redis = get_redis_connection()
    status_key = f"{_STATUS_PREFIX}{job_id}"
    result_key = f"{_RESULT_PREFIX}{job_id}"
    summary_key = f"{_SUMMARY_PREFIX}{job_id}"
    freqs_key = f"{_FREQS_PREFIX}{job_id}"
    lease_key = f"{_LEASE_PREFIX}{job_id}"
    ttl = settings.result_ttl_seconds
    payload_json = json.dumps(result)

    for _ in range(3):
        pipe = redis.pipeline()
        try:
            pipe.watch(status_key, lease_key)
            lease = _get_lease(pipe, job_id)
            if not lease:
                raise PermissionError("lease not found")
            if lease.get("worker_id") != worker_id or lease.get("lease_id") != lease_id:
                raise PermissionError("lease mismatch")

            status = _decode_map(pipe.hgetall(status_key))
            if status.get("status") == "finished":
                existing_result = pipe.get(result_key)
                existing_summary = pipe.get(summary_key)
                existing_freqs = pipe.get(freqs_key)
                pipe.unwatch()
                if (
                    existing_result == payload_json.encode("utf-8")
                    and existing_summary == summary_text.encode("utf-8")
                    and existing_freqs == freqs_csv.encode("utf-8")
                ):
                    return ResultSubmitOutcome(idempotent=True)
                raise ValueError("result already submitted")

            pipe.multi()
            pipe.setex(result_key, ttl, payload_json)
            pipe.setex(summary_key, ttl, summary_text)
            pipe.setex(freqs_key, ttl, freqs_csv)
            pipe.hset(
                status_key,
                mapping={"status": "finished", "detail": "", "updated_at": _now_iso()},
            )
            pipe.expire(status_key, ttl)
            pipe.delete(lease_key)
            pipe.zrem(_LEASE_INDEX, job_id)
            pipe.execute()
            return ResultSubmitOutcome(idempotent=False)
        except WatchError:
            continue
        finally:
            pipe.reset()
    raise RuntimeError("failed to submit result due to concurrent updates")


def submit_failure(
    *,
    job_id: str,
    worker_id: str,
    lease_id: str,
    error_code: str,
    error_message: str,
    traceback: Optional[str] = None,
) -> FailureOutcome:
    settings = get_zpe_settings()
    redis = get_redis_connection()
    status_key = f"{_STATUS_PREFIX}{job_id}"
    lease_key = f"{_LEASE_PREFIX}{job_id}"
    retry_key = f"{_RETRY_PREFIX}{job_id}"
    ttl = settings.result_ttl_seconds

    for _ in range(3):
        pipe = redis.pipeline()
        try:
            pipe.watch(lease_key, retry_key)
            lease = _get_lease(pipe, job_id)
            if not lease:
                raise PermissionError("lease not found")
            if lease.get("worker_id") != worker_id or lease.get("lease_id") != lease_id:
                raise PermissionError("lease mismatch")

            current_retry = pipe.get(retry_key)
            retry_count = int(current_retry or 0) + 1
            detail = f"{error_code}: {error_message}"

            pipe.multi()
            pipe.set(retry_key, retry_count)
            if retry_count <= settings.retry_max:
                delay = min(
                    settings.retry_max_delay_seconds,
                    settings.retry_base_delay_seconds * (2 ** (retry_count - 1)),
                )
                run_at = int(datetime.now(timezone.utc).timestamp()) + int(delay)
                pipe.zadd(_DELAY_ZSET, {job_id: run_at})
                pipe.hset(
                    status_key,
                    mapping={
                        "status": "queued",
                        "detail": f"requeue in {delay}s ({detail})",
                        "updated_at": _now_iso(),
                    },
                )
                pipe.expire(status_key, ttl)
                pipe.expire(retry_key, ttl)
                pipe.delete(lease_key)
                pipe.zrem(_LEASE_INDEX, job_id)
                pipe.execute()
                return FailureOutcome(requeued=True, retry_count=retry_count)

            pipe.lpush(_DLQ, job_id)
            pipe.hset(
                status_key,
                mapping={
                    "status": "failed",
                    "detail": detail,
                    "updated_at": _now_iso(),
                },
            )
            pipe.expire(status_key, ttl)
            pipe.expire(retry_key, ttl)
            pipe.delete(lease_key)
            pipe.zrem(_LEASE_INDEX, job_id)
            pipe.execute()
            return FailureOutcome(requeued=False, retry_count=retry_count)
        except WatchError:
            continue
        finally:
            pipe.reset()
    raise RuntimeError("failed to submit failure due to concurrent updates")
