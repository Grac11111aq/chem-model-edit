from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime, timedelta, timezone
import json
import secrets
from typing import Dict, Optional, cast
from uuid import uuid4

from redis import Redis

from .queue import get_redis_connection
from .settings import get_zpe_settings


_ENROLL_PREFIX = "zpe:enroll:"
_SERVER_PREFIX = "zpe:compute:server:"


def _now() -> datetime:
    return datetime.now(timezone.utc)


def _now_iso() -> str:
    return _now().isoformat()


def _expires_iso(ttl_seconds: int) -> str:
    return (_now() + timedelta(seconds=ttl_seconds)).isoformat()


@dataclass
class EnrollToken:
    token: str
    expires_at: str
    ttl_seconds: int
    label: Optional[str] = None


@dataclass
class ComputeServerRegistration:
    server_id: str
    registered_at: str
    name: Optional[str] = None
    meta: Dict[str, str] = field(default_factory=dict)


class ComputeEnrollStore:
    def __init__(self, redis: Optional[Redis] = None) -> None:
        self.redis = redis or get_redis_connection()

    def create_token(
        self,
        *,
        ttl_seconds: Optional[int] = None,
        label: Optional[str] = None,
    ) -> EnrollToken:
        settings = get_zpe_settings()
        ttl = ttl_seconds if ttl_seconds is not None else settings.enroll_token_ttl_seconds
        if ttl <= 0:
            raise ValueError("ttl_seconds must be >= 1")
        token = secrets.token_urlsafe(32)
        payload = {
            "created_at": _now_iso(),
            "label": label or "",
        }
        key = f"{_ENROLL_PREFIX}{token}"
        self.redis.hset(key, mapping=payload)
        self.redis.expire(key, ttl)
        return EnrollToken(
            token=token,
            expires_at=_expires_iso(ttl),
            ttl_seconds=ttl,
            label=label,
        )

    def consume_token(
        self,
        token: str,
        *,
        name: Optional[str] = None,
        meta: Optional[Dict[str, str]] = None,
    ) -> ComputeServerRegistration:
        key = f"{_ENROLL_PREFIX}{token}"
        deleted = cast(int, self.redis.delete(key))
        if deleted == 0:
            raise KeyError("token not found")
        server_id = f"compute-{uuid4().hex}"
        now = _now_iso()
        payload = {
            "registered_at": now,
            "name": name or "",
        }
        if meta:
            payload["meta"] = json.dumps(meta, ensure_ascii=False)
        self.redis.hset(f"{_SERVER_PREFIX}{server_id}", mapping=payload)
        return ComputeServerRegistration(
            server_id=server_id,
            registered_at=now,
            name=name,
            meta=meta or {},
        )


def get_enroll_store() -> ComputeEnrollStore:
    return ComputeEnrollStore()
