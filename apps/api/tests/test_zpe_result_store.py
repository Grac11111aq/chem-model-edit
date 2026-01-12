from __future__ import annotations

import fakeredis
import pytest

from services.zpe.result_store import RedisResultStore


def test_redis_result_store_roundtrip():
    fake = fakeredis.FakeRedis()
    store = RedisResultStore(redis=fake)

    store.set_status("job-1", "queued")
    status = store.get_status("job-1")
    assert status.status == "queued"
    assert status.detail is None
    assert status.updated_at is not None

    result = {"zpe_ev": 0.1, "freqs_cm": [100.0]}
    store.set_result("job-1", result, summary_text="summary", freqs_csv="freqs")
    assert store.get_result("job-1") == result
    assert store.get_file("job-1", "summary") == "summary"
    assert store.get_file("job-1", "freqs") == "freqs"


def test_redis_result_store_missing_keys():
    fake = fakeredis.FakeRedis()
    store = RedisResultStore(redis=fake)

    with pytest.raises(KeyError):
        store.get_status("missing")

    with pytest.raises(KeyError):
        store.get_result("missing")

    with pytest.raises(KeyError):
        store.get_file("missing", "summary")

    with pytest.raises(ValueError):
        store.get_file("job-1", "unknown")
