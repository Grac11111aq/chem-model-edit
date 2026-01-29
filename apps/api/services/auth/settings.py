from __future__ import annotations

from functools import lru_cache
from typing import Optional

from pydantic_settings import BaseSettings, SettingsConfigDict

from services.zpe.settings import _resolve_env_file, get_zpe_settings


class AuthSettings(BaseSettings):
    model_config = SettingsConfigDict(
        env_prefix="AUTH_",
        env_file=_resolve_env_file(),
        extra="ignore",
    )

    redis_url: Optional[str] = None
    session_ttl_seconds: int = 60 * 60 * 24 * 7
    password_iterations: int = 210_000
    password_pepper: Optional[str] = None


@lru_cache
def get_auth_settings() -> AuthSettings:
    settings = AuthSettings()
    if not settings.redis_url:
        settings.redis_url = get_zpe_settings().redis_url
    return settings
