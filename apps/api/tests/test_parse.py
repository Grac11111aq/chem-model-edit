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
 H 0.0 0.0 0.0
 H 0.0 0.0 0.74
""".strip()


def test_parse_qe_ok():
    response = CLIENT.post('/parse', json={'content': QE_INPUT})
    assert response.status_code == 200
    data = response.json()
    atoms = data['structure']['atoms']
    assert len(atoms) == 2
    assert atoms[0]['symbol'] == 'H'


def test_parse_qe_invalid():
    response = CLIENT.post('/parse', json={'content': 'invalid'})
    assert response.status_code == 400
