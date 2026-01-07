from __future__ import annotations

from fastapi import FastAPI, Request
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse

from models import (
    DeltaTransplantRequest,
    DeltaTransplantResponse,
    ExportRequest,
    ExportResponse,
    LatticeConvertFromParamsRequest,
    LatticeConvertFromVectorsRequest,
    LatticeConvertResponse,
    ParseRequest,
    ParseResponse,
    SupercellRequest,
    SupercellResponse,
    TiledSupercellRequest,
)
from services.export import export_qe_in
from services.lattice import params_to_vectors, vectors_to_params
from services.parse import parse_qe_in
from services.supercell import generate_supercell, generate_tiled_supercell
from services.transplant import transplant_delta

app = FastAPI(title="Chem Model API", version="0.1.0")

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


@app.exception_handler(ValueError)
def handle_value_error(_: Request, exc: ValueError) -> JSONResponse:
    return JSONResponse(status_code=400, content={"detail": str(exc)})


@app.get("/health")
def health() -> dict[str, str]:
    return {"status": "ok"}


@app.post("/parse", response_model=ParseResponse)
def parse_qe(request: ParseRequest) -> ParseResponse:
    structure = parse_qe_in(request.content)
    return ParseResponse(structure=structure)


@app.post("/export", response_model=ExportResponse)
def export_qe(request: ExportRequest) -> ExportResponse:
    content = export_qe_in(request.structure)
    return ExportResponse(content=content)


@app.post("/transplant/delta", response_model=DeltaTransplantResponse)
def transplant_delta_route(
    request: DeltaTransplantRequest,
) -> DeltaTransplantResponse:
    content = transplant_delta(
        request.small_in,
        request.small_out,
        request.large_in,
    )
    return DeltaTransplantResponse(content=content)


@app.post("/supercell", response_model=SupercellResponse)
def supercell(request: SupercellRequest) -> SupercellResponse:
    structure, meta = generate_supercell(
        request.structureA,
        request.structureB,
        request.sequence,
        request.lattice,
    )
    return SupercellResponse(structure=structure, meta=meta)


@app.post("/supercell/tiled", response_model=SupercellResponse)
def supercell_tiled(request: TiledSupercellRequest) -> SupercellResponse:
    structure, meta = generate_tiled_supercell(
        request.structureA,
        request.structureB,
        request.pattern,
        request.lattice,
        check_overlap=request.checkOverlap,
        overlap_tolerance=request.overlapTolerance,
    )
    return SupercellResponse(structure=structure, meta=meta)


@app.post("/lattice/vectors-to-params", response_model=LatticeConvertResponse)
def lattice_vectors_to_params(
    request: LatticeConvertFromVectorsRequest,
) -> LatticeConvertResponse:
    params = vectors_to_params(request.lattice)
    return LatticeConvertResponse(
        lattice=request.lattice, params=params, unit=request.unit
    )


@app.post("/lattice/params-to-vectors", response_model=LatticeConvertResponse)
def lattice_params_to_vectors(
    request: LatticeConvertFromParamsRequest,
) -> LatticeConvertResponse:
    lattice = params_to_vectors(request.params)
    return LatticeConvertResponse(
        lattice=lattice, params=request.params, unit=request.unit
    )
