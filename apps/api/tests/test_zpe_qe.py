from __future__ import annotations

from services.zpe import qe as qe_service


def test_build_espresso_profile_with_argv(monkeypatch):
    class DummyProfile:
        def __init__(self, command: str, pseudo_dir: str, argv: list[str]):
            self.command = command
            self.pseudo_dir = pseudo_dir
            self.argv = argv

    monkeypatch.setattr(qe_service, "EspressoProfile", DummyProfile)
    monkeypatch.setattr(qe_service, "resolve_pw_command", lambda _settings: "pw.x")

    settings = qe_service.ZPESettings(mpi_cmd="mpirun", pw_command="pw.x")
    profile = qe_service.build_espresso_profile(
        pseudo_dir="pseudo",
        use_mpi=True,
        np_core=4,
        environ=True,
        settings=settings,
    )
    assert profile.command == "mpirun"
    assert profile.argv == ["-np", "4", "pw.x", "-environ"]


def test_build_espresso_profile_command_only(monkeypatch):
    class DummyProfile:
        def __init__(self, command: str, pseudo_dir: str):
            self.command = command
            self.pseudo_dir = pseudo_dir

    monkeypatch.setattr(qe_service, "EspressoProfile", DummyProfile)
    monkeypatch.setattr(qe_service, "resolve_pw_command", lambda _settings: "pw.x")

    settings = qe_service.ZPESettings(mpi_cmd="mpirun", pw_command="pw.x")
    profile = qe_service.build_espresso_profile(
        pseudo_dir="pseudo",
        use_mpi=True,
        np_core=4,
        environ=False,
        settings=settings,
    )
    assert profile.command == "mpirun -np 4 pw.x"
