from __future__ import annotations

from typing import List, Optional

from pydantic import BaseModel, Field


class Atom(BaseModel):
    symbol: str = Field(..., min_length=1, max_length=3)
    x: float
    y: float
    z: float


class Structure(BaseModel):
    atoms: List[Atom]
    lattice: Optional[Lattice] = None


class LatticeParams(BaseModel):
    a: float
    b: float
    c: float
    alpha: float
    beta: float
    gamma: float


class LatticeConvertFromVectorsRequest(BaseModel):
    lattice: Lattice
    unit: str = "angstrom"


class LatticeConvertFromParamsRequest(BaseModel):
    params: LatticeParams
    unit: str = "angstrom"


class LatticeConvertResponse(BaseModel):
    lattice: Lattice
    params: LatticeParams
    unit: str


class ParseRequest(BaseModel):
    content: str
    format: Optional[str] = None


class ParseResponse(BaseModel):
    structure: Structure


class StructureCreateRequest(BaseModel):
    content: str
    format: Optional[str] = None


class StructureCreateResponse(BaseModel):
    structure_id: str
    structure: Structure
    source: str


class ExportRequest(BaseModel):
    structure: Structure
    format: Optional[str] = None


class ExportResponse(BaseModel):
    content: str


class DeltaTransplantRequest(BaseModel):
    small_in: str
    small_out: str
    large_in: str


class DeltaTransplantResponse(BaseModel):
    content: str


class Vector3(BaseModel):
    x: float
    y: float
    z: float


class Lattice(BaseModel):
    a: Vector3
    b: Vector3
    c: Vector3


class SupercellRequest(BaseModel):
    structureA: Structure
    structureB: Structure
    sequence: str
    lattice: Lattice


class TiledSupercellRequest(BaseModel):
    structureA: Structure
    structureB: Structure
    pattern: List[List[str]]
    lattice: Lattice
    checkOverlap: bool = False
    overlapTolerance: Optional[float] = None


class SupercellMeta(BaseModel):
    na: int
    nb: int
    layers: int
    overlapCount: int = 0


class SupercellResponse(BaseModel):
    structure: Structure
    meta: SupercellMeta
