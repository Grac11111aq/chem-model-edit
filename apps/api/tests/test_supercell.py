from __future__ import annotations

from ase import Atoms as ASEAtoms

from services.supercell import (
    build_supercell_from_grid,
    generate_supercell,
    generate_tiled_supercell,
)
from models import Atom, Lattice, Structure, SupercellGrid, SupercellGridAxis, Vector3


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


def test_build_supercell_from_grid_default_axis():
    base_atoms = ASEAtoms(
        symbols=["H"],
        positions=[(0.0, 0.0, 0.0)],
        cell=[(1.0, 0.0, 0.0), (0.0, 2.0, 0.0), (0.0, 0.0, 3.0)],
        pbc=True,
    )
    tile_a = ASEAtoms(
        symbols=["H"],
        positions=[(0.0, 0.0, 0.0)],
        cell=base_atoms.get_cell(),
        pbc=True,
    )
    tile_b = ASEAtoms(
        symbols=["He"],
        positions=[(0.5, 0.0, 0.0)],
        cell=base_atoms.get_cell(),
        pbc=True,
    )

    grid = SupercellGrid(
        rows=2,
        cols=2,
        tiles=[["a", "b"], ["a", "b"]],
    )
    atoms_out, overlap_count = build_supercell_from_grid(
        base_atoms,
        grid,
        {"a": tile_a, "b": tile_b},
    )

    assert overlap_count == 0
    assert atoms_out.get_cell()[0][0] == 2.0
    assert atoms_out.get_cell()[1][1] == 4.0
    positions = atoms_out.get_positions()
    assert positions[0][0] == 0.0
    assert positions[1][0] == 1.5
    assert positions[2][1] == 2.0


def test_build_supercell_from_grid_axis_swap():
    base_atoms = ASEAtoms(
        symbols=["H"],
        positions=[(0.0, 0.0, 0.0)],
        cell=[(1.0, 0.0, 0.0), (0.0, 2.0, 0.0), (0.0, 0.0, 3.0)],
        pbc=True,
    )
    grid = SupercellGrid(
        rows=3,
        cols=2,
        tiles=[["a", "a"], ["a", "a"], ["a", "a"]],
        axis=SupercellGridAxis(row="a", col="b"),
    )
    atoms_out, _ = build_supercell_from_grid(
        base_atoms,
        grid,
        {"a": base_atoms},
    )
    assert atoms_out.get_cell()[0][0] == 3.0
    assert atoms_out.get_cell()[1][1] == 4.0
