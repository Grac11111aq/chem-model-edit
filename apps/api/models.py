from __future__ import annotations

from typing import Any, Dict, List, Literal, Optional, Tuple

from pydantic import BaseModel, Field, model_validator


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


class StructureGetResponse(BaseModel):
    structure: Structure


class StructureRegisterRequest(BaseModel):
    structure: Structure
    source: Optional[str] = None


class StructureRegisterResponse(BaseModel):
    structureId: str
    source: str


class StructureImportRequest(BaseModel):
    content: str
    format: Optional[str] = None


class StructureImportResponse(BaseModel):
    structureId: str
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


class SupercellGridAxis(BaseModel):
    row: Literal["a", "b"]
    col: Literal["a", "b"]

    @model_validator(mode="after")
    def check_distinct_axes(self) -> "SupercellGridAxis":
        if self.row == self.col:
            raise ValueError("row and col must map to different axes")
        return self


class SupercellGrid(BaseModel):
    rows: int = Field(..., ge=1)
    cols: int = Field(..., ge=1)
    tiles: List[List[str]]
    axis: Optional[SupercellGridAxis] = None

    @model_validator(mode="after")
    def check_tile_dimensions(self) -> "SupercellGrid":
        if len(self.tiles) != self.rows:
            raise ValueError(f"tiles has {len(self.tiles)} rows, expected {self.rows}")
        for index, row in enumerate(self.tiles):
            if len(row) != self.cols:
                raise ValueError(
                    f"tiles[{index}] has {len(row)} cols, expected {self.cols}"
                )
        return self


class SupercellBuildOptions(BaseModel):
    checkOverlap: bool = False
    overlapTolerance: Optional[float] = None
    validateLattice: Literal["none", "warn", "error"] = "none"


class SupercellBuildOutput(BaseModel):
    includeStructure: bool = False


class SupercellBuildRequest(BaseModel):
    baseStructureId: str
    grid: SupercellGrid
    options: Optional[SupercellBuildOptions] = None
    output: Optional[SupercellBuildOutput] = None


class SupercellBuildMeta(BaseModel):
    rows: int
    cols: int
    tileCount: int
    overlapCount: Optional[int] = None
    baseStructureId: str
    structureIdsUsed: List[str]


class SupercellBuildResponse(BaseModel):
    structureId: str
    structure: Optional[Structure] = None
    meta: SupercellBuildMeta


class SupercellMeta(BaseModel):
    na: int
    nb: int
    layers: int
    overlapCount: int = 0


class SupercellResponse(BaseModel):
    structure: Structure
    meta: SupercellMeta


class ZPEParseRequest(BaseModel):
    content: str
    structure_id: Optional[str] = None


class ZPEParseResponse(BaseModel):
    structure: Structure
    fixed_indices: List[int]
    atomic_species: Dict[str, str] = Field(default_factory=dict)
    kpoints: Optional[Tuple[int, int, int]] = None


class ZPEJobRequest(BaseModel):
    content: str
    mobile_indices: List[int]
    use_environ: bool = False
    input_dir: Optional[str] = None
    calc_mode: Literal["new", "continue"] = "continue"
    structure_id: Optional[str] = None


class ZPEJobResponse(BaseModel):
    job_id: str


class ZPEJobStatusResponse(BaseModel):
    status: str
    detail: Optional[str] = None
    updated_at: Optional[str] = None


class ZPEResult(BaseModel):
    freqs_cm: List[float]
    zpe_ev: float
    s_vib_jmol_k: float
    mobile_indices: List[int]
    fixed_indices: List[int]
    kpts: Tuple[int, int, int]
    delta: float
    low_cut_cm: float
    temperature: float
    use_environ: bool
    qe_input: str
    pseudo_dir: str
    calc_start_time: str
    calc_end_time: str
    elapsed_seconds: float
    cache_checked: int
    cache_deleted: int
    ecutwfc: Optional[float] = None
    ecutrho: Optional[float] = None


class ZPEJobResultResponse(BaseModel):
    result: ZPEResult


class ZPEEnrollTokenRequest(BaseModel):
    ttl_seconds: Optional[int] = None
    label: Optional[str] = None


class ZPEEnrollTokenResponse(BaseModel):
    token: str
    expires_at: str
    ttl_seconds: int
    label: Optional[str] = None


class ZPEComputeRegisterRequest(BaseModel):
    token: str
    name: Optional[str] = None
    meta: Dict[str, Any] = Field(default_factory=dict)


class ZPEComputeRegisterResponse(BaseModel):
    server_id: str
    registered_at: str
    name: Optional[str] = None
    worker_token: str
    token_expires_at: str
    token_ttl_seconds: int


class ZPEComputeRevokeResponse(BaseModel):
    revoked_count: int


class ZPEComputeLeaseResponse(BaseModel):
    job_id: str
    payload: Dict[str, Any]
    lease_id: str
    lease_ttl_seconds: int


class ZPEComputeResultRequest(BaseModel):
    lease_id: str
    result: Dict[str, Any]
    summary_text: str
    freqs_csv: str
    meta: Dict[str, Any] = Field(default_factory=dict)


class ZPEComputeResultResponse(BaseModel):
    ok: bool = True
    idempotent: bool = False


class ZPEComputeFailedRequest(BaseModel):
    lease_id: str
    error_code: str
    error_message: str
    traceback: Optional[str] = None


class ZPEComputeFailedResponse(BaseModel):
    ok: bool = True
    requeued: bool
    retry_count: int
