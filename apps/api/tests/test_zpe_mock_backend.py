from __future__ import annotations

import fakeredis

from models import ZPEJobRequest
from services.zpe import backends
from services.zpe.result_store import RedisResultStore
from services.zpe.settings import ZPESettings


def test_mock_backend_produces_result(monkeypatch):
    fake = fakeredis.FakeRedis()
    store = RedisResultStore(redis=fake)
    settings = ZPESettings(compute_mode="mock", result_store="redis")

    monkeypatch.setattr(backends, "get_result_store", lambda: store)
    monkeypatch.setattr(backends, "get_zpe_settings", lambda: settings)

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
    assert job_id.startswith("mock-")

    status = store.get_status(job_id)
    assert status.status == "finished"

    result = store.get_result(job_id)
    assert result["zpe_ev"] >= 0.0
    assert len(result["freqs_cm"]) == 3
    assert result["mobile_indices"] == [0]

    summary = store.get_file(job_id, "summary")
    assert "ZPE summary" in summary

    freqs_csv = store.get_file(job_id, "freqs")
    assert "frequency_cm^-1" in freqs_csv
