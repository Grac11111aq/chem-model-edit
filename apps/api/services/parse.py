from __future__ import annotations

from io import StringIO
from typing import List

from ase import Atoms as ASEAtoms
from ase.io import read as ase_read
from pymatgen.io.pwscf import PWInput

from models import Atom, Structure


def _from_ase(content: str) -> Structure:
    atoms_result = ase_read(StringIO(content), format="espresso-in")
    if isinstance(atoms_result, list):
        if not atoms_result:
            raise ValueError("ASEが構造を返しませんでした。")
        atoms_obj: ASEAtoms = atoms_result[0]
    else:
        atoms_obj = atoms_result
    symbols = atoms_obj.get_chemical_symbols()
    positions = atoms_obj.get_positions()
    parsed: List[Atom] = []
    for symbol, (x, y, z) in zip(symbols, positions):
        parsed.append(Atom(symbol=symbol, x=float(x), y=float(y), z=float(z)))
    return Structure(atoms=parsed)


def _from_pymatgen(content: str) -> Structure:
    pw_input = PWInput.from_str(content)
    structure = pw_input.structure
    parsed: List[Atom] = []
    for site in structure.sites:
        parsed.append(
            Atom(
                symbol=str(site.specie),
                x=float(site.coords[0]),
                y=float(site.coords[1]),
                z=float(site.coords[2]),
            )
        )
    return Structure(atoms=parsed)


def parse_qe_in(content: str) -> Structure:
    try:
        return _from_ase(content)
    except Exception as ase_error:  # pragma: no cover - fallback path varies
        try:
            return _from_pymatgen(content)
        except Exception as pmg_error:
            raise ValueError(
                f"QE .in のパースに失敗しました: ASE={ase_error}, pymatgen={pmg_error}"
            ) from pmg_error
