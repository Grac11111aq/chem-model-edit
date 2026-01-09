from __future__ import annotations

from pathlib import Path

import fakeredis
from fastapi.testclient import TestClient
from rq import Queue
from rq.job import Job

import main
from services import zpe as zpe_service
from services import zpe_worker
from services.zpe import qe as zpe_qe

QE_INPUT_KPOINTS = """
&CONTROL
  calculation = 'scf'
  prefix = 'zpe_kpts'
  pseudo_dir = './pseudo'
  outdir = './tmp'
/
&SYSTEM
  ibrav = 0
  nat = 1
  ntyp = 1
  ecutwfc = 10.0
  ecutrho = 40.0
/
&ELECTRONS
  conv_thr = 1.0d-8
/
ATOMIC_SPECIES
  H 1.0079 H.pbe-rrkjus_psl.1.0.0.UPF
CELL_PARAMETERS angstrom
  5.0 0.0 0.0
  0.0 5.0 0.0
  0.0 0.0 5.0
ATOMIC_POSITIONS angstrom
  H 0.0 0.0 0.0 1 1 1
K_POINTS (automatic)
  2 3 4 0 0 0
""".strip()


class DummyVibrations:
    def __init__(self, atoms, indices, delta, name):
        self.atoms = atoms
        self.indices = indices
        self.delta = delta
        self.name = name

    def run(self) -> None:
        return None

    def get_frequencies(self):
        return [100.0, 200.0]


class DummyEspresso:
    def __init__(self, **_kwargs):
        pass


def test_zpe_kpoints_reflected(monkeypatch, tmp_path: Path) -> None:
    fake_redis = fakeredis.FakeRedis()
    queue = Queue("zpe", connection=fake_redis, is_async=False)
    monkeypatch.setattr(main, "get_queue", lambda: queue)
    monkeypatch.setattr(
        main,
        "fetch_job",
        lambda job_id: Job.fetch(job_id, connection=fake_redis),
    )

    monkeypatch.setattr(zpe_worker, "Vibrations", DummyVibrations)
    monkeypatch.setattr(zpe_worker, "Espresso", DummyEspresso)
    monkeypatch.setattr(zpe_qe, "resolve_pw_command", lambda _settings: "pw.x")

    monkeypatch.setenv("ZPE_PSEUDO_DIR", str(tmp_path / "pseudo"))
    monkeypatch.setenv("ZPE_WORK_DIR", str(tmp_path / "jobs"))

    # Refresh cached settings to pick up env vars
    zpe_service.get_zpe_settings.cache_clear()

    client = TestClient(main.app)
    resp = client.post(
        "/calc/zpe/jobs",
        json={
            "content": QE_INPUT_KPOINTS,
            "mobile_indices": [0],
            "low_cut_cm": 1.0e9,
        },
    )
    assert resp.status_code == 200
    job_id = resp.json()["job_id"]

    resp = client.get(f"/calc/zpe/jobs/{job_id}/result")
    assert resp.status_code == 200
    result = resp.json()["result"]
    assert result["kpts"] == [2, 3, 4]

    resp = client.get(f"/calc/zpe/jobs/{job_id}/files", params={"kind": "summary"})
    assert resp.status_code == 200
    assert resp.text.startswith("# ZPE summary")
