from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime, timedelta, timezone
import json
import secrets
from typing import Any, Mapping, Optional, cast
from uuid import uuid4

from redis import Redis, WatchError

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
    owner_id: Optional[str] = None


@dataclass
class ComputeServerRegistration:
    server_id: str
    registered_at: str
    name: Optional[str] = None
    meta: Mapping[str, Any] = field(default_factory=dict)
    owner_id: Optional[str] = None


class ComputeEnrollStore:
    def __init__(self, redis: Optional[Redis] = None) -> None:
        self.redis = redis or get_redis_connection()

    def create_token(
        self,
        *,
        ttl_seconds: Optional[int] = None,
        label: Optional[str] = None,
        owner_id: Optional[str] = None,
    ) -> EnrollToken:
        settings = get_zpe_settings()
        ttl = (
            ttl_seconds
            if ttl_seconds is not None
            else settings.enroll_token_ttl_seconds
        )
        if ttl <= 0:
            raise ValueError("ttl_seconds must be >= 1")
        token = secrets.token_urlsafe(32)
        payload = {
            "created_at": _now_iso(),
            "label": label or "",
            "owner_id": owner_id or "",
        }
        key = f"{_ENROLL_PREFIX}{token}"
        pipe = self.redis.pipeline(transaction=True)
        pipe.hset(key, mapping=payload)
        pipe.expire(key, ttl)
        _, expire_ok = pipe.execute()
        if not expire_ok:
            self.redis.delete(key)
            raise RuntimeError("failed to set enroll token expiry")
        return EnrollToken(
            token=token,
            expires_at=_expires_iso(ttl),
            ttl_seconds=ttl,
            label=label,
            owner_id=owner_id,
        )

    def consume_token(
        self,
        token: str,
        *,
        name: Optional[str] = None,
        meta: Optional[Mapping[str, Any]] = None,
    ) -> ComputeServerRegistration:
        key = f"{_ENROLL_PREFIX}{token}"
        server_id = f"compute-{uuid4().hex}"
        now = _now_iso()
        owner_id = None
        payload = {
            "registered_at": now,
            "name": name or "",
        }
        if meta:
            payload["meta"] = json.dumps(meta, ensure_ascii=False)
        server_key = f"{_SERVER_PREFIX}{server_id}"
        pipe = self.redis.pipeline(transaction=True)
        for _ in range(5):
            try:
                pipe.watch(key)
                if not pipe.exists(key):
                    pipe.reset()
                    raise KeyError("token not found")
                token_payload = cast(dict[bytes, bytes], pipe.hgetall(key))
                if token_payload:
                    raw_owner = token_payload.get(b"owner_id") or b""
                    owner_id = raw_owner.decode("utf-8") or None
                pipe.multi()
                pipe.delete(key)
                pipe.hset(server_key, mapping=payload)
                pipe.execute()
                break
            except WatchError:
                pipe.reset()
                continue
        else:
            raise RuntimeError("failed to consume enroll token")
        return ComputeServerRegistration(
            server_id=server_id,
            registered_at=now,
            name=name,
            meta=meta or {},
            owner_id=owner_id,
        )


def get_enroll_store() -> ComputeEnrollStore:
    return ComputeEnrollStore()
