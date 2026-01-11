from __future__ import annotations

import logging
from typing import Literal

from fastapi import FastAPI, HTTPException, Query, Request
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse, Response
from pydantic import ValidationError
from redis.exceptions import RedisError

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
    ZPEJobRequest,
    ZPEJobResponse,
    ZPEJobResultResponse,
    ZPEJobStatusResponse,
    ZPEParseRequest,
    ZPEParseResponse,
    ZPEResult,
)
from services.export import export_qe_in
from services.lattice import params_to_vectors, vectors_to_params
from services.parse import parse_qe_in
from services.supercell import generate_supercell, generate_tiled_supercell
from services.transplant import transplant_delta
from services.zpe import ensure_mobile_indices, enqueue_zpe_job, get_result_store
from services.zpe.parse import parse_atomic_species, parse_kpoints_automatic, parse_qe_structure

app = FastAPI(title="Chem Model API", version="0.1.0")
logger = logging.getLogger(__name__)

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


@app.exception_handler(RedisError)
def handle_redis_error(_: Request, exc: RedisError) -> JSONResponse:
    logger.error("Redis error", exc_info=exc)
    return JSONResponse(
        status_code=503,
        content={"detail": "redis unavailable"},
    )


def _ensure_job_finished(job_id: str) -> None:
    store = get_result_store()
    try:
        status = store.get_status(job_id)
    except KeyError as exc:
        raise HTTPException(status_code=404, detail="job not found") from exc
    if status.status == "failed":
        detail = status.detail or "unknown"
        raise HTTPException(status_code=409, detail=f"job failed: {detail}")
    if status.status != "finished":
        raise HTTPException(
            status_code=409,
            detail=f"job not finished (status={status.status})",
        )


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


@app.post("/calc/zpe/parse", response_model=ZPEParseResponse)
def zpe_parse(request: ZPEParseRequest) -> ZPEParseResponse:
    structure, fixed_indices = parse_qe_structure(request.content)
    kpts = parse_kpoints_automatic(request.content)
    return ZPEParseResponse(
        structure=structure,
        fixed_indices=fixed_indices,
        atomic_species=parse_atomic_species(request.content),
        kpoints=kpts[0] if kpts else None,
    )


@app.post("/calc/zpe/jobs", response_model=ZPEJobResponse)
def zpe_jobs(request: ZPEJobRequest) -> ZPEJobResponse:
    structure, fixed_indices = parse_qe_structure(request.content)
    mobile_indices = ensure_mobile_indices(
        request.mobile_indices, len(structure.atoms), fixed_indices
    )
    payload = request.model_dump()
    payload["mobile_indices"] = mobile_indices
    job_id = enqueue_zpe_job(payload)
    return ZPEJobResponse(job_id=job_id)

@app.get("/calc/zpe/jobs/{job_id}", response_model=ZPEJobStatusResponse)
def zpe_job_status(job_id: str) -> ZPEJobStatusResponse:
    store = get_result_store()
    try:
        status = store.get_status(job_id)
    except KeyError as exc:
        raise HTTPException(status_code=404, detail="job not found") from exc
    return ZPEJobStatusResponse(
        status=status.status,
        detail=status.detail,
        updated_at=status.updated_at,
    )


@app.get("/calc/zpe/jobs/{job_id}/result", response_model=ZPEJobResultResponse)
def zpe_job_result(job_id: str) -> ZPEJobResultResponse:
    store = get_result_store()
    _ensure_job_finished(job_id)
    try:
        result_dict = store.get_result(job_id)
    except KeyError as exc:
        raise HTTPException(status_code=500, detail="result missing after completion") from exc
    try:
        result = ZPEResult(**result_dict)
    except ValidationError as exc:
        raise HTTPException(status_code=500, detail="result data invalid") from exc
    return ZPEJobResultResponse(result=result)


@app.get("/calc/zpe/jobs/{job_id}/files")
def zpe_job_files(
    job_id: str, kind: Literal["summary", "freqs"] = Query(...)
) -> Response:
    store = get_result_store()
    _ensure_job_finished(job_id)
    try:
        payload = store.get_file(job_id, kind)
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc
    except KeyError as exc:
        raise HTTPException(status_code=500, detail="file missing after completion") from exc
    filename = "summary.txt" if kind == "summary" else "freqs.csv"
    media_type = "text/plain" if kind == "summary" else "text/csv"
    return Response(
        content=payload,
        media_type=media_type,
        headers={"Content-Disposition": f'attachment; filename=\"{filename}\"'},
    )
