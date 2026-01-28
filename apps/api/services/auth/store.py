from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timedelta, timezone
import json
import secrets
from typing import Optional, cast
from uuid import uuid4

from redis import Redis, WatchError

from .passwords import hash_password, verify_password
from .settings import get_auth_settings


_USER_PREFIX = "auth:user:"
_EMAIL_PREFIX = "auth:user_email:"
_SESSION_PREFIX = "auth:session:"


def _now() -> datetime:
    return datetime.now(timezone.utc)


def _now_iso() -> str:
    return _now().isoformat()


def _expires_iso(ttl_seconds: int) -> str:
    return (_now() + timedelta(seconds=ttl_seconds)).isoformat()


@dataclass
class AuthUser:
    user_id: str
    email: str
    password_hash: str
    created_at: str


@dataclass
class AuthSession:
    token: str
    user_id: str
    expires_at: str


class AuthStore:
    def __init__(self, redis: Optional[Redis] = None) -> None:
        if redis is not None:
            self.redis = redis
            return
        settings = get_auth_settings()
        if settings.redis_url is None:
            raise RuntimeError("AUTH_REDIS_URL must be set")
        self.redis = Redis.from_url(settings.redis_url)

    def create_user(self, email: str, password: str) -> AuthUser:
        normalized = email.strip().lower()
        if not normalized:
            raise ValueError("email required")
        if len(password) < 8:
            raise ValueError("password must be at least 8 characters")

        settings = get_auth_settings()
        password_hash = hash_password(
            password,
            iterations=settings.password_iterations,
            pepper=settings.password_pepper,
        )
        user_id = f"user-{uuid4().hex}"
        user = AuthUser(
            user_id=user_id,
            email=normalized,
            password_hash=password_hash,
            created_at=_now_iso(),
        )

        email_key = f"{_EMAIL_PREFIX}{normalized}"
        user_key = f"{_USER_PREFIX}{user_id}"
        payload = json.dumps(user.__dict__)

        pipe = self.redis.pipeline(transaction=True)
        for _ in range(5):
            try:
                pipe.watch(email_key)
                if pipe.exists(email_key):
                    pipe.reset()
                    raise ValueError("email already registered")
                pipe.multi()
                pipe.set(email_key, user_id)
                pipe.set(user_key, payload)
                pipe.execute()
                break
            except WatchError:
                pipe.reset()
                continue
        else:
            raise RuntimeError("failed to create user")

        return user

    def get_user_by_id(self, user_id: str) -> Optional[AuthUser]:
        raw = cast(Optional[bytes], self.redis.get(f"{_USER_PREFIX}{user_id}"))
        if not raw:
            return None
        return AuthUser(**json.loads(raw))

    def get_user_by_email(self, email: str) -> Optional[AuthUser]:
        normalized = email.strip().lower()
        user_id = cast(Optional[bytes], self.redis.get(f"{_EMAIL_PREFIX}{normalized}"))
        if not user_id:
            return None
        return self.get_user_by_id(user_id.decode("utf-8"))

    def authenticate(self, email: str, password: str) -> Optional[AuthUser]:
        user = self.get_user_by_email(email)
        if not user:
            return None
        settings = get_auth_settings()
        if not verify_password(
            password, user.password_hash, pepper=settings.password_pepper
        ):
            return None
        return user

    def create_session(self, user_id: str) -> AuthSession:
        settings = get_auth_settings()
        token = secrets.token_urlsafe(32)
        ttl = settings.session_ttl_seconds
        key = f"{_SESSION_PREFIX}{token}"
        ok = self.redis.setex(key, ttl, user_id)
        if not ok:
            raise RuntimeError("failed to create session")
        return AuthSession(token=token, user_id=user_id, expires_at=_expires_iso(ttl))

    def get_user_by_session(
        self, token: str, *, refresh: bool = True
    ) -> Optional[AuthSession]:
        if not token:
            return None
        settings = get_auth_settings()
        key = f"{_SESSION_PREFIX}{token}"
        user_id = cast(Optional[bytes], self.redis.get(key))
        if not user_id:
            return None
        if refresh:
            self.redis.expire(key, settings.session_ttl_seconds)
        return AuthSession(
            token=token,
            user_id=user_id.decode("utf-8"),
            expires_at=_expires_iso(settings.session_ttl_seconds),
        )

    def delete_session(self, token: str) -> None:
        if not token:
            return
        self.redis.delete(f"{_SESSION_PREFIX}{token}")


def get_auth_store() -> AuthStore:
    return AuthStore()
