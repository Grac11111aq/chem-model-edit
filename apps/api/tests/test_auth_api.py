from __future__ import annotations

import fakeredis
from fastapi.testclient import TestClient

import main
from services.auth.store import AuthStore


def _patch_auth(monkeypatch):
    fake = fakeredis.FakeRedis()
    store = AuthStore(fake)
    monkeypatch.setattr(main, "get_auth_store", lambda: store)
    return store


def test_auth_register_login_logout(monkeypatch):
    _patch_auth(monkeypatch)
    client = TestClient(main.app)

    response = client.post(
        "/auth/register",
        json={"email": "user@example.com", "password": "password123"},
    )
    assert response.status_code == 200
    token = response.json()["token"]

    me = client.get("/auth/me", headers={"Authorization": f"Bearer {token}"})
    assert me.status_code == 200
    assert me.json()["user"]["email"] == "user@example.com"

    logout = client.post("/auth/logout", headers={"Authorization": f"Bearer {token}"})
    assert logout.status_code == 200

    after = client.get("/auth/me", headers={"Authorization": f"Bearer {token}"})
    assert after.status_code == 401


def test_auth_login_invalid(monkeypatch):
    _patch_auth(monkeypatch)
    client = TestClient(main.app)

    response = client.post(
        "/auth/register",
        json={"email": "user@example.com", "password": "password123"},
    )
    assert response.status_code == 200

    bad = client.post(
        "/auth/login",
        json={"email": "user@example.com", "password": "wrong"},
    )
    assert bad.status_code == 401


def test_auth_register_duplicate(monkeypatch):
    _patch_auth(monkeypatch)
    client = TestClient(main.app)

    response = client.post(
        "/auth/register",
        json={"email": "user@example.com", "password": "password123"},
    )
    assert response.status_code == 200

    duplicate = client.post(
        "/auth/register",
        json={"email": "user@example.com", "password": "password123"},
    )
    assert duplicate.status_code == 400
