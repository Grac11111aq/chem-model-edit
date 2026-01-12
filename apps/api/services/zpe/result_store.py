from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass
from datetime import datetime, timezone
import json
from typing import Any, Dict, Optional, cast

from redis import Redis

from .queue import get_redis_connection
from .settings import get_zpe_settings


@dataclass
class ZPEJobStatus:
    status: str
    detail: Optional[str] = None
    updated_at: Optional[str] = None


_STATUS_PREFIX = "zpe:status:"
_RESULT_PREFIX = "zpe:result:"
_SUMMARY_PREFIX = "zpe:summary:"
_FREQS_PREFIX = "zpe:freqs:"


def _now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def _decode(value: Optional[bytes]) -> Optional[str]:
    if value is None:
        return None
    return value.decode("utf-8")


def _decode_dict(data: Dict[bytes, bytes]) -> Dict[str, str]:
    return {k.decode("utf-8"): v.decode("utf-8") for k, v in data.items()}


class ResultStore(ABC):
    @abstractmethod
    def set_status(self, job_id: str, status: str, detail: Optional[str] = None) -> None:
        ...

    @abstractmethod
    def get_status(self, job_id: str) -> ZPEJobStatus:
        ...

    @abstractmethod
    def set_result(
        self,
        job_id: str,
        result: Dict[str, Any],
        *,
        summary_text: str,
        freqs_csv: str,
    ) -> None:
        ...

    @abstractmethod
    def get_result(self, job_id: str) -> Dict[str, Any]:
        ...

    @abstractmethod
    def get_file(self, job_id: str, kind: str) -> str:
        ...


class RedisResultStore(ResultStore):
    def __init__(self, redis: Optional[Redis] = None) -> None:
        self.redis = redis or get_redis_connection()
        self.ttl = get_zpe_settings().result_ttl_seconds

    def set_status(self, job_id: str, status: str, detail: Optional[str] = None) -> None:
        key = f"{_STATUS_PREFIX}{job_id}"
        payload = {
            "status": status,
            "detail": detail or "",
            "updated_at": _now_iso(),
        }
        self.redis.hset(key, mapping=payload)
        self.redis.expire(key, self.ttl)

    def get_status(self, job_id: str) -> ZPEJobStatus:
        key = f"{_STATUS_PREFIX}{job_id}"
        data = cast(Dict[bytes, bytes], self.redis.hgetall(key))
        if not data:
            raise KeyError("status not found")
        decoded = _decode_dict(data)
        return ZPEJobStatus(
            status=decoded.get("status", "unknown"),
            detail=decoded.get("detail") or None,
            updated_at=decoded.get("updated_at"),
        )

    def set_result(
        self,
        job_id: str,
        result: Dict[str, Any],
        *,
        summary_text: str,
        freqs_csv: str,
    ) -> None:
        self.redis.setex(f"{_RESULT_PREFIX}{job_id}", self.ttl, json.dumps(result))
        self.redis.setex(f"{_SUMMARY_PREFIX}{job_id}", self.ttl, summary_text)
        self.redis.setex(f"{_FREQS_PREFIX}{job_id}", self.ttl, freqs_csv)

    def get_result(self, job_id: str) -> Dict[str, Any]:
        raw = cast(Optional[bytes], self.redis.get(f"{_RESULT_PREFIX}{job_id}"))
        if raw is None:
            raise KeyError("result not found")
        return json.loads(raw)

    def get_file(self, job_id: str, kind: str) -> str:
        if kind == "summary":
            key = f"{_SUMMARY_PREFIX}{job_id}"
        elif kind == "freqs":
            key = f"{_FREQS_PREFIX}{job_id}"
        else:
            raise ValueError("kind must be 'summary' or 'freqs'.")
        raw = cast(Optional[bytes], self.redis.get(key))
        if raw is None:
            raise KeyError("file not found")
        return _decode(raw) or ""


def get_result_store() -> ResultStore:
    settings = get_zpe_settings()
    if settings.result_store != "redis":
        raise ValueError("result_store must be 'redis'.")
    return RedisResultStore()
