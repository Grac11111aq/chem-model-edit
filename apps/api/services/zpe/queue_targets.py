from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
import json
from typing import Optional, cast
from uuid import uuid4

from redis import Redis

from .queue import get_redis_connection


_TARGET_PREFIX = "zpe:queue:target:"
_USER_TARGETS_PREFIX = "zpe:queue:targets:"
_USER_ACTIVE_PREFIX = "zpe:queue:active:"


def _now() -> datetime:
    return datetime.now(timezone.utc)


def _now_iso() -> str:
    return _now().isoformat()


@dataclass
class QueueTarget:
    target_id: str
    user_id: str
    queue_name: str
    server_id: str
    registered_at: str
    name: Optional[str] = None


class QueueTargetStore:
    def __init__(self, redis: Optional[Redis] = None) -> None:
        self.redis = redis or get_redis_connection()

    def add_target(
        self,
        *,
        user_id: str,
        queue_name: str,
        server_id: str,
        name: Optional[str] = None,
    ) -> QueueTarget:
        target_id = f"qt-{uuid4().hex}"
        target = QueueTarget(
            target_id=target_id,
            user_id=user_id,
            queue_name=queue_name,
            server_id=server_id,
            registered_at=_now_iso(),
            name=name,
        )
        key = f"{_TARGET_PREFIX}{target_id}"
        list_key = f"{_USER_TARGETS_PREFIX}{user_id}"
        payload = json.dumps(target.__dict__)
        pipe = self.redis.pipeline(transaction=True)
        pipe.set(key, payload)
        pipe.rpush(list_key, target_id)
        pipe.execute()
        return target

    def list_targets(self, user_id: str) -> list[QueueTarget]:
        list_key = f"{_USER_TARGETS_PREFIX}{user_id}"
        ids = cast(list[bytes], self.redis.lrange(list_key, 0, -1))
        targets: list[QueueTarget] = []
        for raw in ids:
            target = self.get_target(raw.decode("utf-8"))
            if target:
                targets.append(target)
        return targets

    def get_target(self, target_id: str) -> Optional[QueueTarget]:
        raw = cast(Optional[bytes], self.redis.get(f"{_TARGET_PREFIX}{target_id}"))
        if not raw:
            return None
        return QueueTarget(**json.loads(raw))

    def set_active_target(self, user_id: str, target_id: str) -> None:
        key = f"{_USER_ACTIVE_PREFIX}{user_id}"
        self.redis.set(key, target_id)

    def get_active_target(self, user_id: str) -> Optional[QueueTarget]:
        key = f"{_USER_ACTIVE_PREFIX}{user_id}"
        raw = cast(Optional[bytes], self.redis.get(key))
        if not raw:
            return None
        return self.get_target(raw.decode("utf-8"))

    def ensure_target_owner(self, user_id: str, target_id: str) -> QueueTarget:
        target = self.get_target(target_id)
        if not target or target.user_id != user_id:
            raise KeyError("target not found")
        return target


def get_queue_target_store() -> QueueTargetStore:
    return QueueTargetStore()
