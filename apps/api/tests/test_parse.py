from __future__ import annotations

import pytest
from fastapi.testclient import TestClient

from main import app
from services import parse as parse_service

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
    response = CLIENT.post("/parse", json={"content": QE_INPUT})
    assert response.status_code == 200
    data = response.json()
    atoms = data["structure"]["atoms"]
    assert len(atoms) == 2
    assert atoms[0]["symbol"] == "H"
    lattice = data["structure"].get("lattice")
    assert lattice is not None
    assert lattice["a"]["x"] == 5.0
    assert lattice["b"]["y"] == 5.0
    assert lattice["c"]["z"] == 5.0


def test_parse_qe_invalid():
    response = CLIENT.post("/parse", json={"content": "invalid"})
    assert response.status_code == 400


def test_parse_qe_non_structure_message():
    response = CLIENT.post("/parse", json={"content": "&INPUTPP\\n/\\n"})
    assert response.status_code == 400
    detail = response.json().get("detail", "")
    assert "構造データではありません" in detail


def test_parse_qe_no_manual_fallback(monkeypatch):
    qe_input = """
&CONTROL
  calculation='scf'
/
&SYSTEM
  ibrav=2, celldm(1)=10.2, nat=2, ntyp=1
/
&ELECTRONS
/
ATOMIC_SPECIES
 C 12.0107 C.pbe-rrkjus.UPF
ATOMIC_POSITIONS angstrom
 C 0.0 0.0 0.0
 C 1.0 1.0 1.0
""".strip()

    def _raise(*_args, **_kwargs):
        raise ValueError("force fallback")

    monkeypatch.setattr(parse_service, "_from_ase", _raise)
    monkeypatch.setattr(parse_service, "_from_pymatgen", _raise)

    with pytest.raises(ValueError) as excinfo:
        parse_service.parse_qe_in(qe_input)

    message = str(excinfo.value)
    assert "QE .in のパースに失敗しました" in message
    assert "ASE=force fallback" in message
    assert "pymatgen=force fallback" in message
