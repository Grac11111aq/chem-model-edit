from __future__ import annotations

from typing import Iterable, List, Tuple

from ase import Atoms as ASEAtoms

from models import (
    Atom,
    Lattice,
    Structure,
    SupercellGrid,
    SupercellGridAxis,
    SupercellMeta,
    Vector3,
)


def _apply_shift(
    atoms: Iterable[Atom], shift: Tuple[float, float, float]
) -> List[Atom]:
    dx, dy, dz = shift
    return [
        Atom(
            symbol=atom.symbol,
            x=atom.x + dx,
            y=atom.y + dy,
            z=atom.z + dz,
        )
        for atom in atoms
    ]


def _parse_sequence(sequence: str) -> List[List[str]]:
    blocks = [block.strip() for block in sequence.upper().split(",") if block.strip()]
    return [list(block) for block in blocks]


def _scale_vector(vec: Vector3, factor: int) -> Vector3:
    return Vector3(x=vec.x * factor, y=vec.y * factor, z=vec.z * factor)


def _vector_to_tuple(vec: Vector3) -> Tuple[float, float, float]:
    return (vec.x, vec.y, vec.z)


def generate_supercell(
    structure_a: Structure,
    structure_b: Structure,
    sequence: str,
    lattice: Lattice,
) -> tuple[Structure, SupercellMeta]:
    blocks = _parse_sequence(sequence)
    if not blocks:
        raise ValueError("シーケンスが空です。")

    na = max((len(block) for block in blocks), default=1)
    nb = len(blocks)
    atoms_out: List[Atom] = []
    va = _vector_to_tuple(lattice.a)
    vb = _vector_to_tuple(lattice.b)
    b_shift = (0.0, 0.0, 0.0)
    for block in blocks:
        a_shift = (0.0, 0.0, 0.0)
        for layer in block:
            if layer not in {"A", "B"}:
                raise ValueError(f"未知のシーケンス記号: {layer}")
            source = structure_a if layer == "A" else structure_b
            atoms_out.extend(
                _apply_shift(
                    source.atoms,
                    (
                        a_shift[0] + b_shift[0],
                        a_shift[1] + b_shift[1],
                        a_shift[2] + b_shift[2],
                    ),
                ),
            )
            a_shift = (
                a_shift[0] + va[0],
                a_shift[1] + va[1],
                a_shift[2] + va[2],
            )
        b_shift = (
            b_shift[0] + vb[0],
            b_shift[1] + vb[1],
            b_shift[2] + vb[2],
        )

    out_lattice = Lattice(
        a=_scale_vector(lattice.a, na),
        b=_scale_vector(lattice.b, nb),
        c=lattice.c,
    )
    meta = SupercellMeta(na=na, nb=nb, layers=sum(len(block) for block in blocks))
    return Structure(atoms=atoms_out, lattice=out_lattice), meta


def generate_tiled_supercell(
    structure_a: Structure,
    structure_b: Structure,
    pattern: List[List[str]],
    lattice: Lattice,
    check_overlap: bool = False,
    overlap_tolerance: float | None = None,
) -> tuple[Structure, SupercellMeta]:
    if not pattern or not pattern[0]:
        raise ValueError("タイルパターンが空です。")

    width = len(pattern[0])
    for row in pattern:
        if len(row) != width:
            raise ValueError("タイルパターンは矩形である必要があります。")

    if check_overlap and (overlap_tolerance is None or overlap_tolerance <= 0):
        raise ValueError("重複チェックの許容誤差が無効です。")

    va = _vector_to_tuple(lattice.a)
    vb = _vector_to_tuple(lattice.b)

    atoms_out: List[Atom] = []
    overlap_count = 0
    tolerance_sq = (overlap_tolerance or 0.0) ** 2

    for row_index, row in enumerate(pattern):
        for col_index, cell in enumerate(row):
            if cell not in {"A", "B"}:
                raise ValueError(f"未知のタイル記号: {cell}")
            source = structure_a if cell == "A" else structure_b
            shift = (
                va[0] * col_index + vb[0] * row_index,
                va[1] * col_index + vb[1] * row_index,
                va[2] * col_index + vb[2] * row_index,
            )
            shifted_atoms = _apply_shift(source.atoms, shift)
            if check_overlap and overlap_tolerance is not None:
                for new_atom in shifted_atoms:
                    for existing in atoms_out:
                        dx = existing.x - new_atom.x
                        dy = existing.y - new_atom.y
                        dz = existing.z - new_atom.z
                        if dx * dx + dy * dy + dz * dz <= tolerance_sq:
                            overlap_count += 1
                            break
            atoms_out.extend(shifted_atoms)

    out_lattice = Lattice(
        a=_scale_vector(lattice.a, width),
        b=_scale_vector(lattice.b, len(pattern)),
        c=lattice.c,
    )
    meta = SupercellMeta(
        na=width,
        nb=len(pattern),
        layers=len(pattern) * width,
        overlapCount=overlap_count,
    )
    return Structure(atoms=atoms_out, lattice=out_lattice), meta


def build_supercell_from_grid(
    base_atoms: ASEAtoms,
    grid: SupercellGrid,
    tiles: dict[str, ASEAtoms],
    *,
    check_overlap: bool = False,
    overlap_tolerance: float | None = None,
) -> tuple[ASEAtoms, int]:
    if check_overlap and (overlap_tolerance is None or overlap_tolerance <= 0):
        raise ValueError("重複チェックの許容誤差が無効です。")

    axis = grid.axis or SupercellGridAxis(row="b", col="a")
    base_cell = base_atoms.get_cell()
    base_a = base_cell[0]
    base_b = base_cell[1]

    row_vec = base_a if axis.row == "a" else base_b
    col_vec = base_a if axis.col == "a" else base_b

    symbols: List[str] = []
    positions: List[Tuple[float, float, float]] = []
    overlap_count = 0
    tolerance_sq = (overlap_tolerance or 0.0) ** 2

    for row_index, row in enumerate(grid.tiles):
        for col_index, structure_id in enumerate(row):
            atoms = tiles[structure_id]
            shift = row_vec * row_index + col_vec * col_index
            for symbol, pos in zip(
                atoms.get_chemical_symbols(), atoms.get_positions(), strict=True
            ):
                new_pos = (
                    float(pos[0] + shift[0]),
                    float(pos[1] + shift[1]),
                    float(pos[2] + shift[2]),
                )
                if check_overlap:
                    for existing in positions:
                        dx = existing[0] - new_pos[0]
                        dy = existing[1] - new_pos[1]
                        dz = existing[2] - new_pos[2]
                        if dx * dx + dy * dy + dz * dz <= tolerance_sq:
                            overlap_count += 1
                            break
                positions.append(new_pos)
                symbols.append(symbol)

    # axis.row and axis.col are guaranteed distinct by SupercellGridAxis validation.
    repeat_a = grid.rows if axis.row == "a" else grid.cols
    repeat_b = grid.rows if axis.row == "b" else grid.cols

    out_cell = base_cell.copy()
    out_cell[0] = base_a * repeat_a
    out_cell[1] = base_b * repeat_b

    atoms_out = ASEAtoms(
        symbols=symbols,
        positions=positions,
        cell=out_cell,
        pbc=base_atoms.get_pbc(),
    )
    return atoms_out, overlap_count
