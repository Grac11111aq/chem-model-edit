from __future__ import annotations

from fastapi.testclient import TestClient

from main import app

CLIENT = TestClient(app)

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


def test_zpe_parse_ok():
    response = CLIENT.post("/calc/zpe/parse", json={"content": QE_INPUT})
    assert response.status_code == 200
    data = response.json()
    atoms = data["structure"]["atoms"]
    fixed = data["fixed_indices"]
    assert len(atoms) == 2
    assert fixed == [0]


def test_zpe_parse_invalid():
    response = CLIENT.post("/calc/zpe/parse", json={"content": "invalid"})
    assert response.status_code == 400
