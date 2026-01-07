from __future__ import annotations

from fastapi.testclient import TestClient

from main import app

CLIENT = TestClient(app)


def test_export_qe_ok():
    payload = {
        "structure": {
            "atoms": [
                {"symbol": "C", "x": 0.1, "y": 0.2, "z": 0.3},
                {"symbol": "H", "x": 1.1, "y": 1.2, "z": 1.3},
            ]
        }
    }
    response = CLIENT.post("/export", json=payload)
    assert response.status_code == 200
    content = response.json()["content"]
    assert "ATOMIC_POSITIONS" in content
    assert "C 0.1000000000 0.2000000000 0.3000000000" in content


def test_export_qe_empty():
    response = CLIENT.post("/export", json={"structure": {"atoms": []}})
    assert response.status_code == 400
