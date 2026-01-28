from __future__ import annotations

import fakeredis
from fastapi.testclient import TestClient

import main
from services.auth.store import AuthStore
from services.zpe import backends as zpe_backends
from services.zpe import enroll as zpe_enroll
from services.zpe import job_owner as zpe_job_owner
from services.zpe import queue as zpe_queue
from services.zpe import queue_targets as zpe_queue_targets
from services.zpe import result_store as zpe_store
from services.zpe import worker_auth as zpe_worker_auth
from services.zpe.result_store import RedisResultStore
from services.zpe.settings import ZPESettings

QE_INPUT = """
&CONTROL
  calculation='scf'
/
&SYSTEM
  ibrav=0, nat=1, ntyp=1
/
&ELECTRONS
/
ATOMIC_SPECIES
 H 1.0079 H.pbe-rrkjus.UPF
CELL_PARAMETERS angstrom
  5.0 0.0 0.0
  0.0 5.0 0.0
  0.0 0.0 5.0
ATOMIC_POSITIONS angstrom
 H 0.0 0.0 0.0 1 1 1
""".strip()


def _patch_redis(monkeypatch):
    fake = fakeredis.FakeRedis()
    monkeypatch.setattr(zpe_store, "get_redis_connection", lambda: fake)
    monkeypatch.setattr(zpe_queue, "get_redis_connection", lambda: fake)
    monkeypatch.setattr(zpe_enroll, "get_redis_connection", lambda: fake)
    monkeypatch.setattr(zpe_worker_auth, "get_redis_connection", lambda: fake)
    monkeypatch.setattr(zpe_queue_targets, "get_redis_connection", lambda: fake)
    monkeypatch.setattr(zpe_job_owner, "get_redis_connection", lambda: fake)
    return fake


def _setup_user_and_target(
    client: TestClient, monkeypatch, fake
) -> tuple[dict[str, str], str]:
    store = AuthStore(fake)
    monkeypatch.setattr(main, "get_auth_store", lambda: store)

    response = client.post(
        "/auth/register",
        json={"email": "user@example.com", "password": "password123"},
    )
    payload = response.json()
    token = payload["token"]
    user_id = payload["user"]["user_id"]
    headers = {"Authorization": f"Bearer {token}"}

    enroll = client.post(
        "/calc/zpe/compute/enroll-tokens",
        json={"ttl_seconds": 60},
        headers=headers,
    )
    enroll_token = enroll.json()["token"]

    client.post(
        "/calc/zpe/compute/servers/register",
        json={
            "token": enroll_token,
            "name": "server-1",
            "queue_name": "zpe",
        },
    )
    return headers, user_id


def test_zpe_mock_api_flow(monkeypatch):
    fake = _patch_redis(monkeypatch)
    store = RedisResultStore(redis=fake)
    settings = ZPESettings(compute_mode="mock", result_store="redis")

    monkeypatch.setattr(zpe_backends, "get_result_store", lambda: store)
    monkeypatch.setattr(zpe_backends, "get_zpe_settings", lambda: settings)

    client = TestClient(main.app)
    headers, _user_id = _setup_user_and_target(client, monkeypatch, fake)
    response = client.post(
        "/calc/zpe/jobs",
        json={
            "content": QE_INPUT,
            "mobile_indices": [0],
            "use_environ": False,
            "input_dir": None,
            "calc_mode": "continue",
        },
        headers=headers,
    )
    assert response.status_code == 200
    job_id = response.json()["job_id"]

    status = client.get(f"/calc/zpe/jobs/{job_id}", headers=headers)
    assert status.status_code == 200
    assert status.json()["status"] == "finished"

    result = client.get(f"/calc/zpe/jobs/{job_id}/result", headers=headers)
    assert result.status_code == 200
    payload = result.json()["result"]
    assert payload["zpe_ev"] >= 0.0
    assert payload["mobile_indices"] == [0]

    summary = client.get(
        f"/calc/zpe/jobs/{job_id}/files",
        params={"kind": "summary"},
        headers=headers,
    )
    assert summary.status_code == 200
    assert "ZPE summary" in summary.text

    freqs = client.get(
        f"/calc/zpe/jobs/{job_id}/files",
        params={"kind": "freqs"},
        headers=headers,
    )
    assert freqs.status_code == 200
    assert freqs.text.startswith("frequency_cm^-1")


def test_zpe_job_status_missing(monkeypatch):
    fake = _patch_redis(monkeypatch)
    client = TestClient(main.app)
    headers, _user_id = _setup_user_and_target(client, monkeypatch, fake)
    response = client.get("/calc/zpe/jobs/missing-job", headers=headers)
    assert response.status_code == 404


def test_zpe_job_result_not_finished(monkeypatch):
    fake = _patch_redis(monkeypatch)
    store = RedisResultStore(redis=fake)
    job_id = "job-not-finished"
    store.set_status(job_id, "queued")

    client = TestClient(main.app)
    headers, user_id = _setup_user_and_target(client, monkeypatch, fake)
    owner_store = zpe_job_owner.JobOwnerStore(redis=fake)
    owner_store.set_owner(job_id, user_id)
    response = client.get(f"/calc/zpe/jobs/{job_id}/result", headers=headers)
    assert response.status_code == 409
