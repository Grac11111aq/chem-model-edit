from __future__ import annotations

from io import StringIO
import re
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from ase import Atoms as ASEAtoms
from ase.io import read as ase_read

from models import Atom, Lattice, Structure, Vector3


_POSITION_HEADER = re.compile(r"^\s*ATOMIC_POSITIONS\b", re.IGNORECASE)
_POSITION_ANY = re.compile(r"^\s*ATOMIC_POSITIONS\b", re.IGNORECASE | re.MULTILINE)
_BLOCK_START = re.compile(
    r"^\s*(CELL_PARAMETERS|K_POINTS|ATOMIC_SPECIES|CONSTRAINTS|OCCUPATIONS|&|/)",
    re.IGNORECASE,
)
_SPECIES_HEADER = re.compile(r"^\s*ATOMIC_SPECIES\b", re.IGNORECASE)


def _iter_atomic_position_lines(content: str) -> List[str]:
    lines = content.splitlines()
    start = None
    for idx, line in enumerate(lines):
        if _POSITION_HEADER.search(line):
            start = idx + 1
            break
    if start is None:
        raise ValueError("ATOMIC_POSITIONS block not found in QE input.")

    atom_lines: List[str] = []
    for raw in lines[start:]:
        stripped = raw.strip()
        if not stripped:
            break
        if _BLOCK_START.match(stripped):
            break
        if stripped.startswith("!"):
            continue
        atom_lines.append(raw)
    return atom_lines


def _is_zero_flag(flag: str) -> bool:
    try:
        return abs(float(flag)) < 1.0e-12
    except ValueError:
        return False


def extract_fixed_indices(content: str) -> List[int]:
    atom_lines = _iter_atomic_position_lines(content)
    fixed: List[int] = []
    atom_index = 0
    for raw in atom_lines:
        line = raw.split("!", 1)[0].strip()
        if not line:
            continue
        tokens = line.split()
        if len(tokens) < 4:
            continue
        flags = tokens[4:7]
        if len(flags) >= 3 and all(_is_zero_flag(flag) for flag in flags[:3]):
            fixed.append(atom_index)
        atom_index += 1
    return fixed


def parse_atomic_species(content: str) -> Dict[str, str]:
    lines = content.splitlines()
    start = None
    for idx, line in enumerate(lines):
        if _SPECIES_HEADER.search(line):
            start = idx + 1
            break
    if start is None:
        return {}
    species: Dict[str, str] = {}
    for raw in lines[start:]:
        stripped = raw.strip()
        if not stripped:
            break
        if _BLOCK_START.match(stripped):
            break
        if stripped.startswith("!"):
            continue
        tokens = stripped.split()
        if len(tokens) < 3:
            continue
        symbol = tokens[0]
        pseudo = tokens[2]
        species[symbol] = pseudo
    return species


def parse_kpoints_automatic(
    content: str,
) -> Optional[Tuple[Tuple[int, int, int], Tuple[int, int, int]]]:
    match = re.search(
        r"K_POINTS\s*\(\s*automatic\s*\)\s*\n\s*([\-\d]+)\s+([\-\d]+)\s+([\-\d]+)\s+([\-\d]+)\s+([\-\d]+)\s+([\-\d]+)",
        content,
        re.IGNORECASE,
    )
    if not match:
        return None
    nkx, nky, nkz, sx, sy, sz = (int(match.group(i)) for i in range(1, 7))
    return (nkx, nky, nkz), (sx, sy, sz)


def _parse_namelist_value(value: str) -> object:
    lowered = value.lower()
    if lowered in (".true.", "true"):
        return True
    if lowered in (".false.", "false"):
        return False
    try:
        if any(ch in lowered for ch in [".", "e", "d"]):
            num = float(lowered.replace("d", "e"))
            if abs(num - int(num)) < 1.0e-12:
                return int(num)
            return num
        return int(lowered)
    except ValueError:
        return value.strip().strip("\"'")


def parse_namelist(content: str, name: str) -> Dict[str, object]:
    match = re.search(rf"&\s*{name}\b(.*?)\n\s*/", content, re.IGNORECASE | re.DOTALL)
    if not match:
        return {}
    body = match.group(1)
    out: Dict[str, object] = {}
    for line in body.splitlines():
        line = line.split("!", 1)[0]
        if not line.strip():
            continue
        for chunk in line.split(","):
            if "=" not in chunk:
                continue
            key, value = chunk.split("=", 1)
            k = key.strip().lower()
            v = value.strip()
            out[k] = _parse_namelist_value(v)
    return out


def parse_qe_atoms(content: str) -> ASEAtoms:
    if not _POSITION_ANY.search(content):
        raise ValueError("構造データではありません。")
    try:
        atoms_result = ase_read(StringIO(content), format="espresso-in")
    except Exception as exc:  # pragma: no cover - ASE内部例外は多様
        raise ValueError(f"QE .in のパースに失敗しました: {exc}") from exc
    if isinstance(atoms_result, list):
        if not atoms_result:
            raise ValueError("ASEが構造を返しませんでした。")
        return atoms_result[0]
    return atoms_result


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


def parse_qe_structure(content: str) -> Tuple[Structure, List[int]]:
    atoms = parse_qe_atoms(content)
    symbols = atoms.get_chemical_symbols()
    positions = atoms.get_positions()
    parsed_atoms: List[Atom] = []
    for symbol, (x, y, z) in zip(symbols, positions):
        parsed_atoms.append(Atom(symbol=symbol, x=float(x), y=float(y), z=float(z)))
    lattice = _lattice_from_vectors(atoms.get_cell().tolist())
    fixed_indices = extract_fixed_indices(content)
    return Structure(atoms=parsed_atoms, lattice=lattice), fixed_indices


def ensure_mobile_indices(
    mobile_indices: Iterable[int],
    natoms: int,
    fixed_indices: Iterable[int],
) -> List[int]:
    mobile = sorted({int(i) for i in mobile_indices})
    if not mobile:
        raise ValueError("可動原子が空です。")
    if any(idx < 0 or idx >= natoms for idx in mobile):
        raise ValueError("可動原子のインデックスが範囲外です。")
    fixed_set = {int(i) for i in fixed_indices}
    if fixed_set.intersection(mobile):
        raise ValueError("固定原子が可動原子に含まれています。")
    return mobile
