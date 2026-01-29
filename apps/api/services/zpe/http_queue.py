from __future__ import annotations

import json
from typing import Any, Dict
from uuid import uuid4

from .queue import get_redis_connection
from .result_store import get_result_store
from .settings import get_zpe_settings


_QUEUE_KEY = "zpe:queue"
_PAYLOAD_PREFIX = "zpe:payload:"


def enqueue_http_job(payload: Dict[str, Any]) -> str:
    settings = get_zpe_settings()
    job_id = f"http-{uuid4().hex}"
    redis = get_redis_connection()
    store = get_result_store()
    payload_key = f"{_PAYLOAD_PREFIX}{job_id}"

    redis.setex(payload_key, settings.result_ttl_seconds, json.dumps(payload))
    redis.lpush(_QUEUE_KEY, job_id)
    store.set_status(job_id, "queued")
    return job_id
