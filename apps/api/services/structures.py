from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timedelta, timezone
from threading import RLock
from typing import Dict, Optional
from uuid import uuid4

from ase import Atoms as ASEAtoms

from models import Structure
from services.cif import atoms_to_cif, cif_to_bcif
from services.parse import parse_qe_atoms, structure_from_ase

DEFAULT_TTL_SECONDS = 3600


@dataclass
class StoredStructure:
    atoms: ASEAtoms
    source: str
    cif: str
    created_at: datetime
    expires_at: datetime


class StructureStore:
    def __init__(self, ttl_seconds: int = DEFAULT_TTL_SECONDS) -> None:
        self._ttl_seconds = ttl_seconds
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
            expires_at=now + timedelta(seconds=self._ttl_seconds),
        )
        with self._lock:
            self._purge_expired(now)
            self._items[structure_id] = entry
        return structure_id

    def get(self, structure_id: str) -> Optional[StoredStructure]:
        now = datetime.now(timezone.utc)
        with self._lock:
            self._purge_expired(now)
            entry = self._items.get(structure_id)
            if not entry:
                return None
            if entry.expires_at <= now:
                self._items.pop(structure_id, None)
                return None
            return entry

    def _purge_expired(self, now: datetime) -> None:
        expired = [key for key, entry in self._items.items() if entry.expires_at <= now]
        for key in expired:
            self._items.pop(key, None)


_STORE = StructureStore()


def create_structure_from_qe(content: str) -> tuple[str, Structure, str]:
    atoms, source = parse_qe_atoms(content)
    cif = atoms_to_cif(atoms)
    structure = structure_from_ase(atoms)
    structure_id = _STORE.create(atoms=atoms, source=source, cif=cif)
    return structure_id, structure, source


def get_structure_bcif(
    structure_id: str,
    *,
    lossy: bool,
    precision: Optional[int],
) -> bytes:
    entry = _STORE.get(structure_id)
    if not entry:
        raise KeyError(structure_id)
    return cif_to_bcif(entry.cif, lossy=lossy, precision=precision)


def get_structure_entry(structure_id: str) -> StoredStructure:
    entry = _STORE.get(structure_id)
    if not entry:
        raise KeyError(structure_id)
    return entry


def register_structure_atoms(atoms: ASEAtoms, source: str) -> str:
    cif = atoms_to_cif(atoms)
    return _STORE.create(atoms=atoms, source=source, cif=cif)
