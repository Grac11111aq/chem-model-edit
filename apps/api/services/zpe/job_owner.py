from __future__ import annotations

from typing import Optional, cast

from redis import Redis

from .queue import get_redis_connection
from .settings import get_zpe_settings


_OWNER_PREFIX = "zpe:job:owner:"


class JobOwnerStore:
    def __init__(self, redis: Optional[Redis] = None) -> None:
        self.redis = redis or get_redis_connection()

    def set_owner(self, job_id: str, user_id: str) -> None:
        settings = get_zpe_settings()
        key = f"{_OWNER_PREFIX}{job_id}"
        ttl = settings.result_ttl_seconds
        ok = self.redis.setex(key, ttl, user_id)
        if not ok:
            raise RuntimeError("failed to set job owner")

    def get_owner(self, job_id: str) -> Optional[str]:
        raw = cast(Optional[bytes], self.redis.get(f"{_OWNER_PREFIX}{job_id}"))
        if not raw:
            return None
        return raw.decode("utf-8")


def get_job_owner_store() -> JobOwnerStore:
    return JobOwnerStore()
