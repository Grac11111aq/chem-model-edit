from __future__ import annotations

import fakeredis
import pytest

from services.zpe.enroll import ComputeEnrollStore


def test_enroll_token_roundtrip():
    fake = fakeredis.FakeRedis()
    store = ComputeEnrollStore(redis=fake)

    token = store.create_token(ttl_seconds=60, label="lab")
    assert token.token
    assert token.ttl_seconds == 60
    assert token.expires_at

    registration = store.consume_token(token.token, name="server-1", meta={"host": "node"})
    assert registration.server_id.startswith("compute-")
    assert registration.name == "server-1"
    assert registration.meta["host"] == "node"

    with pytest.raises(KeyError):
        store.consume_token(token.token)
