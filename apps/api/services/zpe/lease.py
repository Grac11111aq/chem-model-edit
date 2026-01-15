from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timedelta, timezone
import json
from typing import Any, Dict, Optional
from uuid import uuid4

from redis import Redis

from .queue import get_redis_connection
from .result_store import get_result_store
from .settings import get_zpe_settings


_QUEUE_KEY = "zpe:queue"
_PAYLOAD_PREFIX = "zpe:payload:"
_LEASE_PREFIX = "zpe:lease:"
_LEASE_INDEX = "zpe:lease:index"
_DELAY_ZSET = "zpe:delay"


def _now() -> datetime:
    return datetime.now(timezone.utc)


def _now_iso() -> str:
    return _now().isoformat()


@dataclass
class Lease:
    job_id: str
    payload: Dict[str, Any]
    lease_id: str
    lease_ttl_seconds: int


def lease_next_job(worker_id: str) -> Optional[Lease]:
    settings = get_zpe_settings()
    redis = get_redis_connection()
    store = get_result_store()

    _reap_expired_leases(redis, store)
    _promote_due_jobs(redis)

    lease_id = uuid4().hex
    expires_at = _now() + timedelta(seconds=settings.lease_ttl_seconds)
    expires_iso = expires_at.isoformat()
    expires_ts = int(expires_at.timestamp())
    lease_ttl = int(settings.lease_ttl_seconds)

    script = (
        "local job_id = redis.call('RPOP', KEYS[1])\n"
        "if not job_id then return nil end\n"
        "local lease_key = ARGV[5] .. job_id\n"
        "redis.call('HSET', lease_key, 'worker_id', ARGV[1], 'lease_id', ARGV[2], 'expires_at', ARGV[3])\n"
        "redis.call('EXPIRE', lease_key, ARGV[4])\n"
        "redis.call('ZADD', KEYS[2], ARGV[5], job_id)\n"
        "return job_id"
    )
    job_id_raw = redis.eval(
        script,
        2,
        _QUEUE_KEY,
        _LEASE_INDEX,
        worker_id,
        lease_id,
        expires_iso,
        lease_ttl,
        expires_ts,
        _LEASE_PREFIX,
    )
    if not job_id_raw:
        return None

    job_id = job_id_raw.decode("utf-8") if isinstance(job_id_raw, bytes) else str(job_id_raw)
    payload_raw = redis.get(f"{_PAYLOAD_PREFIX}{job_id}")
    if payload_raw is None:
        redis.delete(f"{_LEASE_PREFIX}{job_id}")
        redis.zrem(_LEASE_INDEX, job_id)
        store.set_status(job_id, "failed", detail="payload missing")
        return None

    try:
        payload = json.loads(payload_raw)
    except json.JSONDecodeError:
        redis.delete(f"{_LEASE_PREFIX}{job_id}")
        redis.zrem(_LEASE_INDEX, job_id)
        store.set_status(job_id, "failed", detail="payload corrupted")
        return None
    return Lease(
        job_id=job_id,
        payload=payload,
        lease_id=lease_id,
        lease_ttl_seconds=settings.lease_ttl_seconds,
    )


def _reap_expired_leases(redis: Redis, store: Any) -> None:
    now_ts = int(_now().timestamp())
    expired = redis.zrangebyscore(_LEASE_INDEX, 0, now_ts)
    if not expired:
        return
    for job_id_raw in expired:
        job_id = job_id_raw.decode("utf-8") if isinstance(job_id_raw, bytes) else str(job_id_raw)
        lease_key = f"{_LEASE_PREFIX}{job_id}"
        lease = redis.hgetall(lease_key)
        if not lease:
            redis.zrem(_LEASE_INDEX, job_id)
            continue
        redis.delete(lease_key)
        redis.zrem(_LEASE_INDEX, job_id)
        redis.lpush(_QUEUE_KEY, job_id)
        store.set_status(job_id, "queued", detail="lease expired")


def _promote_due_jobs(redis: Redis) -> None:
    now_ts = int(_now().timestamp())
    due = redis.zrangebyscore(_DELAY_ZSET, 0, now_ts)
    if not due:
        return
    pipe = redis.pipeline(transaction=True)
    for job_id in due:
        pipe.lpush(_QUEUE_KEY, job_id)
    pipe.zrem(_DELAY_ZSET, *due)
    pipe.execute()


def get_lease_store(redis: Optional[Redis] = None) -> Redis:
    return redis or get_redis_connection()
