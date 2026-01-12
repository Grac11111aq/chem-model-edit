from __future__ import annotations

from services.supercell import generate_supercell, generate_tiled_supercell
from models import Atom, Lattice, Structure, Vector3


def _make_lattice(
    a: tuple[float, float, float], b: tuple[float, float, float]
) -> Lattice:
    return Lattice(
        a=Vector3(x=a[0], y=a[1], z=a[2]),
        b=Vector3(x=b[0], y=b[1], z=b[2]),
        c=Vector3(x=0.0, y=0.0, z=1.0),
    )


def test_generate_supercell_sequence_shifts():
    lattice = _make_lattice((1.0, 0.0, 0.0), (0.0, 1.0, 0.0))
    structure_a = Structure(atoms=[Atom(symbol="A", x=0.0, y=0.0, z=0.0)])
    structure_b = Structure(atoms=[Atom(symbol="B", x=0.0, y=0.0, z=0.0)])

    structure, meta = generate_supercell(structure_a, structure_b, "AB", lattice)

    assert meta.na == 2
    assert meta.nb == 1
    assert meta.layers == 2
    assert len(structure.atoms) == 2
    assert structure.atoms[0].symbol == "A"
    assert structure.atoms[0].x == 0.0
    assert structure.atoms[1].symbol == "B"
    assert structure.atoms[1].x == 1.0


def test_generate_tiled_supercell_overlap_count():
    lattice = _make_lattice((0.4, 0.0, 0.0), (0.0, 1.0, 0.0))
    structure_a = Structure(atoms=[Atom(symbol="A", x=0.0, y=0.0, z=0.0)])
    structure_b = Structure(atoms=[Atom(symbol="B", x=0.0, y=0.0, z=0.0)])

    structure, meta = generate_tiled_supercell(
        structure_a,
        structure_b,
        [["A", "A"]],
        lattice,
        check_overlap=True,
        overlap_tolerance=0.5,
    )

    assert len(structure.atoms) == 2
    assert meta.na == 2
    assert meta.nb == 1
    assert meta.overlapCount == 1
