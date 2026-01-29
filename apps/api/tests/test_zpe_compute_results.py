from __future__ import annotations

import fakeredis

from services.zpe import compute_results as zpe_results
from services.zpe.settings import ZPESettings


def _patch(monkeypatch):
    fake = fakeredis.FakeRedis()
    monkeypatch.setattr(zpe_results, "get_redis_connection", lambda: fake)
    settings = ZPESettings(result_ttl_seconds=60, retry_max=2, retry_base_delay_seconds=1, retry_max_delay_seconds=5)
    monkeypatch.setattr(zpe_results, "get_zpe_settings", lambda: settings)
    return fake


def test_submit_result_idempotent(monkeypatch):
    fake = _patch(monkeypatch)
    job_id = "job-1"
    lease_id = "lease-1"
    fake.hset(
        f"zpe:lease:{job_id}",
        mapping={"worker_id": "worker-1", "lease_id": lease_id, "expires_at": "x"},
    )

    outcome = zpe_results.submit_result(
        job_id=job_id,
        worker_id="worker-1",
        lease_id=lease_id,
        result={"zpe_ev": 0.1},
        summary_text="summary",
        freqs_csv="freqs",
    )
    assert outcome.idempotent is False

    outcome_again = zpe_results.submit_result(
        job_id=job_id,
        worker_id="worker-1",
        lease_id=lease_id,
        result={"zpe_ev": 0.1},
        summary_text="summary",
        freqs_csv="freqs",
    )
    assert outcome_again.idempotent is True


def test_submit_failure_requeue_and_dlq(monkeypatch):
    fake = _patch(monkeypatch)
    job_id = "job-2"
    lease_id = "lease-2"
    fake.hset(
        f"zpe:lease:{job_id}",
        mapping={"worker_id": "worker-1", "lease_id": lease_id, "expires_at": "x"},
    )

    outcome = zpe_results.submit_failure(
        job_id=job_id,
        worker_id="worker-1",
        lease_id=lease_id,
        error_code="ERR",
        error_message="boom",
    )
    assert outcome.requeued is True
    assert outcome.retry_count == 1
    assert fake.zcard("zpe:delay") == 1

    fake.hset(
        f"zpe:lease:{job_id}",
        mapping={"worker_id": "worker-1", "lease_id": lease_id, "expires_at": "x"},
    )
    outcome2 = zpe_results.submit_failure(
        job_id=job_id,
        worker_id="worker-1",
        lease_id=lease_id,
        error_code="ERR",
        error_message="boom",
    )
    assert outcome2.requeued is True
    assert outcome2.retry_count == 2

    fake.hset(
        f"zpe:lease:{job_id}",
        mapping={"worker_id": "worker-1", "lease_id": lease_id, "expires_at": "x"},
    )
    outcome3 = zpe_results.submit_failure(
        job_id=job_id,
        worker_id="worker-1",
        lease_id=lease_id,
        error_code="ERR",
        error_message="boom",
    )
    assert outcome3.requeued is False
    assert outcome3.retry_count == 3
    assert fake.lrange("zpe:dlq", 0, -1) == [job_id.encode("utf-8")]
