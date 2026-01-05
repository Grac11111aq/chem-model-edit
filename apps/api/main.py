from __future__ import annotations

from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware

from models import (
    ExportRequest,
    ExportResponse,
    ParseRequest,
    ParseResponse,
    SupercellRequest,
    SupercellResponse,
)
from services.export import export_qe_in
from services.parse import parse_qe_in
from services.supercell import generate_supercell

app = FastAPI(title="Chem Model API", version="0.1.0")

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


@app.get("/health")
def health() -> dict[str, str]:
    return {"status": "ok"}


@app.post("/parse", response_model=ParseResponse)
def parse_qe(request: ParseRequest) -> ParseResponse:
    try:
        structure = parse_qe_in(request.content)
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc
    return ParseResponse(structure=structure)


@app.post("/export", response_model=ExportResponse)
def export_qe(request: ExportRequest) -> ExportResponse:
    try:
        content = export_qe_in(request.structure)
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc
    return ExportResponse(content=content)


@app.post("/supercell", response_model=SupercellResponse)
def supercell(request: SupercellRequest) -> SupercellResponse:
    try:
        structure, meta = generate_supercell(
            request.structureA,
            request.structureB,
            request.sequence,
            (request.lattice.a.x, request.lattice.a.y, request.lattice.a.z),
            (request.lattice.b.x, request.lattice.b.y, request.lattice.b.z),
        )
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc
    return SupercellResponse(structure=structure, meta=meta)
