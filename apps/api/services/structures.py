from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
from threading import RLock
from typing import Dict, Optional
from uuid import uuid4

from ase import Atoms as ASEAtoms

from models import Structure
from services.cif import atoms_to_cif
from services.parse import parse_qe_atoms, structure_from_ase


@dataclass
class StoredStructure:
    atoms: ASEAtoms
    source: str
    cif: str
    created_at: datetime


class StructureStore:
    def __init__(self) -> None:
        self._items: Dict[str, StoredStructure] = {}
        self._lock = RLock()

    def create(self, atoms: ASEAtoms, source: str, cif: str) -> str:
        now = datetime.now(timezone.utc)
        structure_id = uuid4().hex
        entry = StoredStructure(
            atoms=atoms,
            source=source,
            cif=cif,
            created_at=now,
        )
        with self._lock:
            self._items[structure_id] = entry
        return structure_id

    def get(self, structure_id: str) -> Optional[StoredStructure]:
        with self._lock:
            entry = self._items.get(structure_id)
            if not entry:
                return None
            return entry


_STORE = StructureStore()


def create_structure_from_qe(content: str) -> tuple[str, Structure, str]:
    atoms, source = parse_qe_atoms(content)
    cif = atoms_to_cif(atoms)
    structure = structure_from_ase(atoms)
    structure_id = _STORE.create(atoms=atoms, source=source, cif=cif)
    return structure_id, structure, source


def get_structure_cif(structure_id: str) -> str:
    entry = _STORE.get(structure_id)
    if not entry:
        raise KeyError(structure_id)
    return entry.cif


def get_structure_entry(structure_id: str) -> StoredStructure:
    entry = _STORE.get(structure_id)
    if not entry:
        raise KeyError(structure_id)
    return entry


def get_structure(structure_id: str) -> Structure:
    entry = get_structure_entry(structure_id)
    return structure_from_ase(entry.atoms)


def register_structure_atoms(atoms: ASEAtoms, source: str) -> str:
    cif = atoms_to_cif(atoms)
    return _STORE.create(atoms=atoms, source=source, cif=cif)
