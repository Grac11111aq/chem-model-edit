from __future__ import annotations

import fakeredis

from models import ZPEJobRequest
from services.zpe import backends
from services.zpe import http_queue as zpe_http
from services.zpe import queue as zpe_queue
from services.zpe import result_store as zpe_store
from services.zpe.settings import ZPESettings


def test_remote_http_enqueue(monkeypatch):
    fake = fakeredis.FakeRedis()
    settings = ZPESettings(compute_mode="remote-http", result_store="redis", result_ttl_seconds=60)

    monkeypatch.setattr(backends, "get_zpe_settings", lambda: settings)
    monkeypatch.setattr(backends, "get_result_store", lambda: zpe_store.RedisResultStore(redis=fake))
    monkeypatch.setattr(zpe_queue, "get_redis_connection", lambda: fake)
    monkeypatch.setattr(zpe_http, "get_redis_connection", lambda: fake)
    monkeypatch.setattr(zpe_http, "get_zpe_settings", lambda: settings)
    monkeypatch.setattr(zpe_store, "get_redis_connection", lambda: fake)

    qe_input = """
&control
/
ATOMIC_POSITIONS angstrom
H 0 0 0 1 1 1
""".strip()

    payload = ZPEJobRequest(
        content=qe_input,
        mobile_indices=[0],
        use_environ=False,
        input_dir=None,
        calc_mode="continue",
    ).model_dump()

    job_id = backends.enqueue_zpe_job(payload)
    assert job_id.startswith("http-")

    assert fake.lrange("zpe:queue", 0, -1) == [job_id.encode("utf-8")]
    assert fake.get(f"zpe:payload:{job_id}") is not None

    status = zpe_store.RedisResultStore(redis=fake).get_status(job_id)
    assert status.status == "queued"
