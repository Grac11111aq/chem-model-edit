from __future__ import annotations

from io import BytesIO

from ase import Atoms as ASEAtoms
from ase.io import write as ase_write


def atoms_to_cif(atoms: ASEAtoms, *, wrap: bool = False) -> str:
    buffer = BytesIO()
    ase_write(buffer, atoms, format="cif", wrap=wrap)
    return buffer.getvalue().decode("utf-8")
