from __future__ import annotations

from dataclasses import dataclass
from io import StringIO
import math
import re
from typing import Any, Callable, List, Tuple

ase_read: Callable[..., Any] | None
try:
    from ase.io import read as ase_read
except Exception:
    ase_read = None


@dataclass
class Atom:
    elem: str
    x: float
    y: float
    z: float
    flags: Tuple[int, int, int]


@dataclass
class Structure:
    cell: List[Tuple[float, float, float]]
    atoms: List[Atom]


_FLOAT_RE = re.compile(r"^[+-]?(?:\d+\.?\d*|\.\d+)(?:[eEdD][+-]?\d+)?$")
_SECTION_STOP = (
    "k_points",
    "cell_parameters",
    "atomic_species",
    "constraints",
    "occupations",
)


def _is_float(tok: str) -> bool:
    return bool(_FLOAT_RE.match(tok))


def _parse_float(tok: str) -> float:
    return float(tok.replace("D", "E").replace("d", "e"))


def _read_cell(lines: List[str]) -> List[Tuple[float, float, float]]:
    for i, line in enumerate(lines):
        if line.strip().lower().startswith("cell_parameters"):
            if i + 3 >= len(lines):
                raise ValueError("CELL_PARAMETERS ブロックが不完全です。")
            v1 = [_parse_float(x) for x in lines[i + 1].split()[:3]]
            v2 = [_parse_float(x) for x in lines[i + 2].split()[:3]]
            v3 = [_parse_float(x) for x in lines[i + 3].split()[:3]]
            return [
                (v1[0], v1[1], v1[2]),
                (v2[0], v2[1], v2[2]),
                (v3[0], v3[1], v3[2]),
            ]
    raise ValueError("CELL_PARAMETERS が見つかりません。")


def _read_atoms_with_flags(lines: List[str]) -> List[Atom]:
    start = None
    for i, line in enumerate(lines):
        if line.strip().lower().startswith("atomic_positions"):
            start = i + 1
            break
    if start is None:
        raise ValueError("ATOMIC_POSITIONS が見つかりません。")

    atoms: List[Atom] = []
    for line in lines[start:]:
        s = line.strip()
        if not s:
            break
        if s.lower().startswith(_SECTION_STOP):
            break
        parts = s.split()
        if len(parts) < 4 or not _is_float(parts[1]):
            break
        if len(parts) < 7 or not all(p in ("0", "1") for p in parts[4:7]):
            raise ValueError("ATOMIC_POSITIONS に 0/1 フラグが必要です。")
        elem = parts[0]
        x, y, z = _parse_float(parts[1]), _parse_float(parts[2]), _parse_float(parts[3])
        flags = (int(parts[4]), int(parts[5]), int(parts[6]))
        atoms.append(Atom(elem, x, y, z, flags))
    if not atoms:
        raise ValueError("ATOMIC_POSITIONS から原子を取得できませんでした。")
    return atoms


def read_qe_input(content: str) -> Structure:
    lines = content.splitlines()
    cell = _read_cell(lines)
    atoms = _read_atoms_with_flags(lines)
    return Structure(cell=cell, atoms=atoms)


def _last_positions_from_output(
    content: str, nat_expected: int
) -> List[Tuple[float, float, float]]:
    if ase_read is not None:
        try:
            ase_atoms = ase_read(StringIO(content), format="espresso-out", index=-1)
            if isinstance(ase_atoms, list):
                if not ase_atoms:
                    raise ValueError("ASEが構造を返しませんでした。")
                ase_atoms = ase_atoms[-1]
            positions = ase_atoms.get_positions()
            if len(positions) == nat_expected:
                return [(float(x), float(y), float(z)) for x, y, z in positions]
        except Exception:
            pass

    lines = content.splitlines()
    idxs = [
        i
        for i, ln in enumerate(lines)
        if ln.strip().lower().startswith("atomic_positions")
    ]
    if not idxs:
        raise ValueError("出力に ATOMIC_POSITIONS が見つかりません。")

    for start in reversed(idxs):
        coords: List[Tuple[float, float, float]] = []
        for ln in lines[start + 1 :]:
            s = ln.strip()
            if not s:
                break
            parts = s.split()
            if len(parts) < 4 or not _is_float(parts[1]):
                break
            coords.append(
                (_parse_float(parts[1]), _parse_float(parts[2]), _parse_float(parts[3]))
            )
            if len(coords) == nat_expected:
                break
        if len(coords) == nat_expected:
            return coords

    raise ValueError(
        "出力の ATOMIC_POSITIONS から最終座標を取得できませんでした。"
        " tprnfor=.true. の出力を確認してください。"
    )


def _wrap_delta(dx: float, length: float) -> float:
    if length <= 0:
        return dx
    return dx - round(dx / length) * length


def _assert_supported_cell(cell: List[Tuple[float, float, float]]) -> None:
    (ax, ay, az), (bx, by, bz) = cell[0], cell[1]
    dot = ax * bx + ay * by + az * bz
    if abs(dot) > 1.0e-6:
        raise ValueError("非直交セルは未対応です。a,b が直交していません。")
    if abs(ay) > 1.0e-6 or abs(az) > 1.0e-6 or abs(bx) > 1.0e-6 or abs(bz) > 1.0e-6:
        raise ValueError("a は x 軸、b は y 軸に整列したセルのみ対応しています。")


def _replace_atomic_positions(base: str, new_atoms: List[Atom]) -> str:
    lines = base.splitlines()

    start = None
    for i, ln in enumerate(lines):
        if ln.strip().lower().startswith("atomic_positions"):
            start = i
            break
    if start is None:
        raise ValueError("ATOMIC_POSITIONS が見つかりません。")

    end = start + 1
    while end < len(lines):
        s = lines[end].strip()
        if not s:
            break
        if s.lower().startswith(_SECTION_STOP):
            break
        parts = s.split()
        if len(parts) < 4 or not _is_float(parts[1]):
            break
        end += 1

    new_block = []
    for atom in new_atoms:
        fx, fy, fz = atom.flags
        new_block.append(
            f"  {atom.elem}  {atom.x:.8f}  {atom.y:.8f}  {atom.z:.8f}  {fx} {fy} {fz}"
        )

    lines = lines[: start + 1] + new_block + lines[end:]
    return "\n".join(lines) + "\n"


def transplant_delta(small_in: str, small_out: str, large_in: str) -> str:
    small0 = read_qe_input(small_in)
    large0 = read_qe_input(large_in)

    _assert_supported_cell(small0.cell)

    Lx = math.sqrt(sum(c * c for c in small0.cell[0]))
    Ly = math.sqrt(sum(c * c for c in small0.cell[1]))
    if Lx <= 0 or Ly <= 0:
        raise ValueError("セル長が不正です。")

    small_rel_xyz = _last_positions_from_output(
        small_out, nat_expected=len(small0.atoms)
    )

    small_mov = [i for i, atom in enumerate(small0.atoms) if atom.flags != (0, 0, 0)]
    large_mov = [i for i, atom in enumerate(large0.atoms) if atom.flags != (0, 0, 0)]

    if len(small_mov) == 0 or len(large_mov) == 0:
        raise ValueError("可動原子が見つかりません。")
    if len(small_mov) != len(large_mov):
        raise ValueError(
            f"可動原子数が一致しません: small={len(small_mov)} large={len(large_mov)}"
        )

    deltas: List[Tuple[float, float, float]] = []
    for idx in small_mov:
        base = small0.atoms[idx]
        xr, yr, zr = small_rel_xyz[idx]
        dx = _wrap_delta(xr - base.x, Lx)
        dy = _wrap_delta(yr - base.y, Ly)
        dz = zr - base.z
        deltas.append((dx, dy, dz))

    new_atoms = [Atom(a.elem, a.x, a.y, a.z, a.flags) for a in large0.atoms]
    for j, idx in enumerate(large_mov):
        dx, dy, dz = deltas[j]
        atom = new_atoms[idx]
        x = (atom.x + dx) % Lx
        y = (atom.y + dy) % Ly
        z = atom.z + dz
        new_atoms[idx] = Atom(atom.elem, x, y, z, atom.flags)

    return _replace_atomic_positions(large_in, new_atoms)
