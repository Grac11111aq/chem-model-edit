from __future__ import annotations

import base64
import hashlib
import secrets
from typing import Tuple


def _hash_pbkdf2(password: str, *, salt: bytes, iterations: int) -> bytes:
    return hashlib.pbkdf2_hmac("sha256", password.encode("utf-8"), salt, iterations)


def hash_password(password: str, *, iterations: int, pepper: str | None = None) -> str:
    salt = secrets.token_bytes(16)
    password_value = f"{password}{pepper or ''}"
    digest = _hash_pbkdf2(password_value, salt=salt, iterations=iterations)
    salt_b64 = base64.b64encode(salt).decode("ascii")
    digest_b64 = base64.b64encode(digest).decode("ascii")
    return f"pbkdf2_sha256${iterations}${salt_b64}${digest_b64}"


def _parse_hash(stored: str) -> Tuple[int, bytes, bytes]:
    try:
        algo, iterations_str, salt_b64, digest_b64 = stored.split("$", 3)
    except ValueError as exc:
        raise ValueError("invalid password hash format") from exc
    if algo != "pbkdf2_sha256":
        raise ValueError("unsupported password hash algorithm")
    iterations = int(iterations_str)
    salt = base64.b64decode(salt_b64)
    digest = base64.b64decode(digest_b64)
    return iterations, salt, digest


def verify_password(password: str, stored: str, *, pepper: str | None = None) -> bool:
    iterations, salt, digest = _parse_hash(stored)
    password_value = f"{password}{pepper or ''}"
    candidate = _hash_pbkdf2(password_value, salt=salt, iterations=iterations)
    return secrets.compare_digest(candidate, digest)
