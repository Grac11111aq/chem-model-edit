from __future__ import annotations

from io import StringIO
import math
import re
from typing import List, Optional, Sequence, Tuple

from ase import Atoms as ASEAtoms
from ase.io import read as ase_read
from pymatgen.io.pwscf import PWInput

from models import Atom, Lattice, Structure, Vector3


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
    lattice = _lattice_from_vectors(atoms_obj.get_cell().tolist())
    return Structure(atoms=parsed, lattice=lattice)


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
    lattice = _lattice_from_vectors(structure.lattice.matrix.tolist())
    return Structure(atoms=parsed, lattice=lattice)


BOHR_TO_ANGSTROM = 0.529177210903
_POSITION_HEADER = re.compile(
    r"^\s*ATOMIC_POSITIONS\b(?:\s*[\{\(]?\s*([a-zA-Z_]+)\s*[\}\)]?)?",
    re.IGNORECASE,
)
_POSITION_ANY = re.compile(r"^\s*ATOMIC_POSITIONS\b", re.IGNORECASE | re.MULTILINE)
_CELL_HEADER = re.compile(
    r"^\s*CELL_PARAMETERS\b(?:\s*[\{\(]?\s*([a-zA-Z_]+)\s*[\}\)]?)?",
    re.IGNORECASE,
)
_SYSTEM_BLOCK = re.compile(r"&SYSTEM(.*?)/", re.IGNORECASE | re.DOTALL)
_SECTION_HEADINGS = {
    "ATOMIC_SPECIES",
    "ATOMIC_POSITIONS",
    "CELL_PARAMETERS",
    "K_POINTS",
    "CONSTRAINTS",
    "ATOMIC_FORCES",
    "OCCUPATIONS",
}


def _strip_inline_comment(line: str) -> str:
    for token in ("!", "#"):
        if token in line:
            line = line.split(token, 1)[0]
    return line.strip()


def _parse_float(value: str) -> float:
    return float(value.replace("D", "E").replace("d", "e"))


def _parse_system_params(content: str) -> dict[str, Optional[float]]:
    block_match = _SYSTEM_BLOCK.search(content)
    block = block_match.group(1) if block_match else content

    def _find_value(pattern: str) -> Optional[float]:
        match = re.search(pattern, block, re.IGNORECASE)
        if not match:
            return None
        return _parse_float(match.group(1))

    def _find_int(pattern: str) -> Optional[int]:
        match = re.search(pattern, block, re.IGNORECASE)
        if not match:
            return None
        return int(match.group(1))

    celldm1 = _find_value(r"celldm\s*\(\s*1\s*\)\s*=\s*([0-9eEdD.+-]+)")
    celldm2 = _find_value(r"celldm\s*\(\s*2\s*\)\s*=\s*([0-9eEdD.+-]+)")
    celldm3 = _find_value(r"celldm\s*\(\s*3\s*\)\s*=\s*([0-9eEdD.+-]+)")
    a_param = _find_value(r"\bA\s*=\s*([0-9eEdD.+-]+)")
    b_param = _find_value(r"\bB\s*=\s*([0-9eEdD.+-]+)")
    c_param = _find_value(r"\bC\s*=\s*([0-9eEdD.+-]+)")
    ibrav = _find_int(r"\bibrav\s*=\s*([+-]?\d+)")

    a_angstrom = a_param if a_param is not None else None
    if a_angstrom is None and celldm1 is not None:
        a_angstrom = celldm1 * BOHR_TO_ANGSTROM

    b_angstrom = b_param
    if b_angstrom is None and a_angstrom is not None and celldm2 is not None:
        b_angstrom = a_angstrom * celldm2

    c_angstrom = c_param
    if c_angstrom is None and a_angstrom is not None and celldm3 is not None:
        c_angstrom = a_angstrom * celldm3

    return {
        "ibrav": int(ibrav) if ibrav is not None else None,
        "a": a_angstrom,
        "b": b_angstrom,
        "c": c_angstrom,
    }


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


def _lattice_from_ibrav(
    params: dict[str, Optional[float]],
) -> Optional[List[Tuple[float, float, float]]]:
    ibrav_val = params.get("ibrav")
    if ibrav_val is None:
        return None
    ibrav = int(ibrav_val)
    a = params.get("a")
    b = params.get("b")
    c = params.get("c")
    if a is None:
        return None

    if ibrav == 1:
        return [(a, 0.0, 0.0), (0.0, a, 0.0), (0.0, 0.0, a)]
    if ibrav == 2:
        return [
            (0.0, a / 2.0, a / 2.0),
            (a / 2.0, 0.0, a / 2.0),
            (a / 2.0, a / 2.0, 0.0),
        ]
    if ibrav == 3:
        return [
            (-a / 2.0, a / 2.0, a / 2.0),
            (a / 2.0, -a / 2.0, a / 2.0),
            (a / 2.0, a / 2.0, -a / 2.0),
        ]

    if ibrav == 4:
        if c is None:
            return None
        return [(a, 0.0, 0.0), (-a / 2.0, a * math.sqrt(3) / 2.0, 0.0), (0.0, 0.0, c)]
    if ibrav == 6:
        if c is None:
            return None
        return [(a, 0.0, 0.0), (0.0, a, 0.0), (0.0, 0.0, c)]
    if ibrav == 7:
        if c is None:
            return None
        return [
            (a / 2.0, -a / 2.0, c / 2.0),
            (a / 2.0, a / 2.0, c / 2.0),
            (-a / 2.0, a / 2.0, c / 2.0),
        ]
    if ibrav == 8:
        if b is None or c is None:
            return None
        return [(a, 0.0, 0.0), (0.0, b, 0.0), (0.0, 0.0, c)]
    if ibrav == 9:
        if b is None or c is None:
            return None
        return [(a / 2.0, b / 2.0, 0.0), (-a / 2.0, b / 2.0, 0.0), (0.0, 0.0, c)]
    if ibrav == 10:
        if b is None or c is None:
            return None
        return [
            (0.0, b / 2.0, c / 2.0),
            (a / 2.0, 0.0, c / 2.0),
            (a / 2.0, b / 2.0, 0.0),
        ]
    if ibrav == 11:
        if b is None or c is None:
            return None
        return [
            (-a / 2.0, b / 2.0, c / 2.0),
            (a / 2.0, -b / 2.0, c / 2.0),
            (a / 2.0, b / 2.0, -c / 2.0),
        ]
    return None


def _parse_cell_parameters(
    content: str, params: dict[str, Optional[float]]
) -> Optional[List[Tuple[float, float, float]]]:
    lines = content.splitlines()
    for idx, raw in enumerate(lines):
        match = _CELL_HEADER.match(raw)
        if not match:
            continue
        unit = (match.group(1) or "alat").strip().lower()
        vectors: List[Tuple[float, float, float]] = []
        for next_line in lines[idx + 1 :]:
            line = _strip_inline_comment(next_line)
            if not line:
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            vectors.append(
                (
                    _parse_float(parts[0]),
                    _parse_float(parts[1]),
                    _parse_float(parts[2]),
                )
            )
            if len(vectors) == 3:
                break
        if len(vectors) < 3:
            raise ValueError("CELL_PARAMETERS のベクトルが不足しています。")

        if unit in ("angstrom", "ang"):
            return vectors
        if unit == "bohr":
            return [
                (
                    v[0] * BOHR_TO_ANGSTROM,
                    v[1] * BOHR_TO_ANGSTROM,
                    v[2] * BOHR_TO_ANGSTROM,
                )
                for v in vectors
            ]
        if unit == "alat":
            a = params.get("a")
            if a is None:
                raise ValueError(
                    "CELL_PARAMETERS が alat 指定ですが格子定数が見つかりません。"
                )
            return [(v[0] * a, v[1] * a, v[2] * a) for v in vectors]
        return vectors
    return None


def _extract_atomic_positions(
    content: str,
) -> Tuple[str, List[Tuple[str, Tuple[float, float, float]]]]:
    lines = content.splitlines()
    for idx, raw in enumerate(lines):
        match = _POSITION_HEADER.match(raw)
        if not match:
            continue
        unit = (match.group(1) or "alat").strip().lower()
        entries: List[Tuple[str, Tuple[float, float, float]]] = []
        for next_line in lines[idx + 1 :]:
            line = _strip_inline_comment(next_line)
            if not line:
                if entries:
                    break
                continue
            if line.startswith("&"):
                break
            first = line.split()[0].upper()
            if entries and first in _SECTION_HEADINGS:
                break
            parts = line.split()
            if len(parts) < 4:
                continue
            symbol = parts[0]
            coords = (
                _parse_float(parts[1]),
                _parse_float(parts[2]),
                _parse_float(parts[3]),
            )
            entries.append((symbol, coords))
        if not entries:
            raise ValueError("ATOMIC_POSITIONS の原子座標が見つかりません。")
        return unit, entries
    raise ValueError("ATOMIC_POSITIONS が見つかりません。")


def _from_manual_positions(content: str) -> Structure:
    unit, entries = _extract_atomic_positions(content)
    params = _parse_system_params(content)
    lattice = _parse_cell_parameters(content, params)
    if lattice is None:
        lattice = _lattice_from_ibrav(params)

    parsed: List[Atom] = []
    if unit in ("crystal", "crystal_sg"):
        if lattice is None:
            raise ValueError(
                "ATOMIC_POSITIONS が crystal 指定ですが格子情報が不足しています。"
            )
        for symbol, coords in entries:
            x = (
                coords[0] * lattice[0][0]
                + coords[1] * lattice[1][0]
                + coords[2] * lattice[2][0]
            )
            y = (
                coords[0] * lattice[0][1]
                + coords[1] * lattice[1][1]
                + coords[2] * lattice[2][1]
            )
            z = (
                coords[0] * lattice[0][2]
                + coords[1] * lattice[1][2]
                + coords[2] * lattice[2][2]
            )
            parsed.append(Atom(symbol=symbol, x=float(x), y=float(y), z=float(z)))
        return Structure(atoms=parsed, lattice=_lattice_from_vectors(lattice))

    if unit == "bohr":
        scale = BOHR_TO_ANGSTROM
    elif unit == "alat":
        alat = params.get("a")
        if alat is None:
            raise ValueError(
                "ATOMIC_POSITIONS が alat 指定ですが格子定数が見つかりません。"
            )
        scale = alat
    else:
        scale = 1.0

    for symbol, coords in entries:
        parsed.append(
            Atom(
                symbol=symbol,
                x=float(coords[0] * scale),
                y=float(coords[1] * scale),
                z=float(coords[2] * scale),
            )
        )
    return Structure(atoms=parsed, lattice=_lattice_from_vectors(lattice))


def parse_qe_in(content: str) -> Structure:
    if not _POSITION_ANY.search(content):
        raise ValueError("構造データではありません。")
    try:
        return _from_ase(content)
    except Exception as ase_error:  # pragma: no cover - fallback path varies
        try:
            return _from_pymatgen(content)
        except Exception as pmg_error:
            try:
                return _from_manual_positions(content)
            except Exception as manual_error:
                raise ValueError(
                    "QE .in のパースに失敗しました: "
                    f"ASE={ase_error}, pymatgen={pmg_error}, manual={manual_error}"
                ) from manual_error
