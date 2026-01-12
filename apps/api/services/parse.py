from __future__ import annotations

from io import StringIO
import re
from typing import List, Optional, Sequence, Tuple

from ase import Atoms as ASEAtoms
from ase.io import read as ase_read
from pymatgen.io.pwscf import PWInput

from models import Atom, Lattice, Structure, Vector3


def _ase_atoms_from_content(content: str) -> ASEAtoms:
    atoms_result = ase_read(StringIO(content), format="espresso-in")
    if isinstance(atoms_result, list):
        if not atoms_result:
            raise ValueError("ASEが構造を返しませんでした。")
        return atoms_result[0]
    return atoms_result


def _pymatgen_atoms_from_content(content: str) -> ASEAtoms:
    pw_input = PWInput.from_str(content)
    structure = pw_input.structure
    symbols = [
        getattr(site.specie, "symbol", str(site.specie))
        for site in structure.sites
    ]
    positions = [site.coords for site in structure.sites]
    return ASEAtoms(
        symbols=symbols,
        positions=positions,
        cell=structure.lattice.matrix.tolist(),
        pbc=True,
    )


def structure_from_ase(atoms_obj: ASEAtoms) -> Structure:
    symbols = atoms_obj.get_chemical_symbols()
    positions = atoms_obj.get_positions()
    parsed: List[Atom] = []
    for symbol, (x, y, z) in zip(symbols, positions):
        parsed.append(Atom(symbol=symbol, x=float(x), y=float(y), z=float(z)))
    lattice = _lattice_from_vectors(atoms_obj.get_cell().tolist())
    return Structure(atoms=parsed, lattice=lattice)


def _from_ase(content: str) -> Structure:
    atoms_obj = _ase_atoms_from_content(content)
    return structure_from_ase(atoms_obj)


def _from_pymatgen(content: str) -> Structure:
    atoms_obj = _pymatgen_atoms_from_content(content)
    return structure_from_ase(atoms_obj)


_POSITION_ANY = re.compile(r"^\s*ATOMIC_POSITIONS\b", re.IGNORECASE | re.MULTILINE)


def _lattice_from_vectors(
    vectors: Optional[Sequence[Sequence[float]]],
) -> Optional[Lattice]:
    if not vectors or len(vectors) < 3:
        return None
    flat = [abs(value) for row in vectors[:3] for value in row[:3]]
    if all(value < 1.0e-8 for value in flat):
        return None
    a, b, c = vectors[0], vectors[1], vectors[2]
    return Lattice(
        a=Vector3(x=float(a[0]), y=float(a[1]), z=float(a[2])),
        b=Vector3(x=float(b[0]), y=float(b[1]), z=float(b[2])),
        c=Vector3(x=float(c[0]), y=float(c[1]), z=float(c[2])),
    )


def parse_qe_in(content: str) -> Structure:
    if not _POSITION_ANY.search(content):
        raise ValueError("構造データではありません。")
    try:
        return _from_ase(content)
    except Exception as ase_error:  # pragma: no cover - fallback path varies
        try:
            return _from_pymatgen(content)
        except Exception as pmg_error:
            raise ValueError(
                f"QE .in のパースに失敗しました: ASE={ase_error}, pymatgen={pmg_error}"
            ) from pmg_error


def parse_qe_atoms(content: str) -> Tuple[ASEAtoms, str]:
    if not _POSITION_ANY.search(content):
        raise ValueError("構造データではありません。")
    try:
        return _ase_atoms_from_content(content), "ase"
    except Exception as ase_error:  # pragma: no cover - fallback path varies
        try:
            return _pymatgen_atoms_from_content(content), "pymatgen"
        except Exception as pmg_error:
            raise ValueError(
                f"QE .in のパースに失敗しました: ASE={ase_error}, pymatgen={pmg_error}"
            ) from pmg_error
