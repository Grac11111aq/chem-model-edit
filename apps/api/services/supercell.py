from __future__ import annotations

from typing import Iterable, List, Tuple

from models import Atom, Structure, SupercellMeta


def _shift_atoms(atoms: Iterable[Atom], shift: Tuple[float, float, float]) -> List[Atom]:
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
    blocks = [block.strip() for block in sequence.upper().split(',') if block.strip()]
    return [list(block) for block in blocks]


def generate_supercell(
    structure_a: Structure,
    structure_b: Structure,
    sequence: str,
    va: Tuple[float, float, float],
    vb: Tuple[float, float, float],
) -> tuple[Structure, SupercellMeta]:
    blocks = _parse_sequence(sequence)
    if not blocks:
        raise ValueError('シーケンスが空です。')

    na = max((len(block) for block in blocks), default=1)
    nb = len(blocks)
    atoms_out: List[Atom] = []
    b_shift = (0.0, 0.0, 0.0)
    for block in blocks:
        a_shift = (0.0, 0.0, 0.0)
        for layer in block:
            if layer not in {'A', 'B'}:
                raise ValueError(f"未知のシーケンス記号: {layer}")
            source = structure_a if layer == 'A' else structure_b
            atoms_out.extend(
                _shift_atoms(
                    source.atoms,
                    (
                        a_shift[0] + b_shift[0],
                        a_shift[1] + b_shift[1],
                        a_shift[2] + b_shift[2],
                    ),
                )
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

    meta = SupercellMeta(na=na, nb=nb, layers=sum(len(block) for block in blocks))
    return Structure(atoms=atoms_out), meta
