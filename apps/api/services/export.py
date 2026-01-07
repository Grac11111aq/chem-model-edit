from __future__ import annotations

from typing import Iterable

from models import Atom, Structure


def _atoms_to_qe_positions(atoms: Iterable[Atom]) -> str:
    lines = ["ATOMIC_POSITIONS angstrom"]
    for atom in atoms:
        lines.append(f"{atom.symbol} {atom.x:.10f} {atom.y:.10f} {atom.z:.10f}")
    return "\n".join(lines)


def export_qe_in(structure: Structure) -> str:
    if not structure.atoms:
        raise ValueError("原子が空です。")
    return _atoms_to_qe_positions(structure.atoms)
