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


class ParseRequest(BaseModel):
    content: str
    format: Optional[str] = None


class ParseResponse(BaseModel):
    structure: Structure


class ExportRequest(BaseModel):
    structure: Structure
    format: Optional[str] = None


class ExportResponse(BaseModel):
    content: str


class Vector3(BaseModel):
    x: float
    y: float
    z: float


class SupercellLattice(BaseModel):
    a: Vector3
    b: Vector3
    c: Vector3


class SupercellRequest(BaseModel):
    structureA: Structure
    structureB: Structure
    sequence: str
    lattice: SupercellLattice


class SupercellMeta(BaseModel):
    na: int
    nb: int
    layers: int


class SupercellResponse(BaseModel):
    structure: Structure
    meta: SupercellMeta
