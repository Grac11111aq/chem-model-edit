from __future__ import annotations

import fakeredis
from fastapi.testclient import TestClient
from rq import Queue

import main

CLIENT = TestClient(main.app)

QE_INPUT = """
&CONTROL
  calculation='scf'
/
&SYSTEM
  ibrav=0, nat=2, ntyp=1
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
 H 0.0 0.0 0.0 0 0 0
 H 0.0 0.0 0.74 1 1 1
""".strip()


def test_zpe_job_enqueue(monkeypatch):
    fake_redis = fakeredis.FakeRedis()
    queue = Queue("zpe", connection=fake_redis)
    monkeypatch.setattr(main, "get_queue", lambda: queue)

    response = CLIENT.post(
        "/calc/zpe/jobs",
        json={"content": QE_INPUT, "mobile_indices": [1]},
    )
    assert response.status_code == 200
    job_id = response.json()["job_id"]
    assert job_id
    assert queue.fetch_job(job_id) is not None


def test_zpe_job_empty_mobile(monkeypatch):
    fake_redis = fakeredis.FakeRedis()
    queue = Queue("zpe", connection=fake_redis)
    monkeypatch.setattr(main, "get_queue", lambda: queue)

    response = CLIENT.post(
        "/calc/zpe/jobs",
        json={"content": QE_INPUT, "mobile_indices": []},
    )
    assert response.status_code == 400
