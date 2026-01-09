from __future__ import annotations

import warnings

from fastapi import FastAPI, HTTPException, Query, Request
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse, JSONResponse

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
)
from services.export import export_qe_in
from services.lattice import params_to_vectors, vectors_to_params
from services.parse import parse_qe_in
from services.supercell import generate_supercell, generate_tiled_supercell
from services.transplant import transplant_delta
from services.zpe import (
    ensure_mobile_indices,
    fetch_job,
    get_queue,
    get_zpe_settings,
    parse_qe_structure,
    resolve_job_file,
)
from services.zpe_worker import run_zpe_job

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


@app.post("/calc/zpe/parse", response_model=ZPEParseResponse)
def zpe_parse(request: ZPEParseRequest) -> ZPEParseResponse:
    structure, fixed_indices = parse_qe_structure(request.content)
    return ZPEParseResponse(structure=structure, fixed_indices=fixed_indices)


@app.post("/calc/zpe/jobs", response_model=ZPEJobResponse)
def zpe_enqueue(request: ZPEJobRequest) -> ZPEJobResponse:
    structure, fixed_indices = parse_qe_structure(request.content)
    ensure_mobile_indices(request.mobile_indices, len(structure.atoms), fixed_indices)
    payload = request.model_dump()
    settings = get_zpe_settings()
    job = get_queue().enqueue(
        run_zpe_job,
        payload,
        job_timeout=settings.job_timeout_seconds,
        result_ttl=settings.result_ttl_seconds,
    )
    return ZPEJobResponse(job_id=job.id)


@app.get("/calc/zpe/jobs/{job_id}", response_model=ZPEJobStatusResponse)
def zpe_status(job_id: str) -> ZPEJobStatusResponse:
    try:
        job = fetch_job(job_id)
    except Exception as exc:  # pragma: no cover - depends on redis behavior
        raise HTTPException(status_code=404, detail="ジョブが見つかりません。") from exc
    status = job.get_status()
    error = job.exc_info if status == "failed" else None
    return ZPEJobStatusResponse(
        job_id=job.id,
        status=status,
        enqueued_at=job.enqueued_at.isoformat() if job.enqueued_at else None,
        started_at=job.started_at.isoformat() if job.started_at else None,
        ended_at=job.ended_at.isoformat() if job.ended_at else None,
        error=error,
    )


@app.get("/calc/zpe/jobs/{job_id}/result", response_model=ZPEJobResultResponse)
def zpe_result(job_id: str) -> ZPEJobResultResponse:
    try:
        job = fetch_job(job_id)
    except Exception as exc:  # pragma: no cover - depends on redis behavior
        raise HTTPException(status_code=404, detail="ジョブが見つかりません。") from exc
    status = job.get_status()
    if status == "failed":
        raise ValueError(
            "ジョブが失敗しています。詳細は /calc/zpe/jobs/{id} を確認してください。"
        )
    if status != "finished":
        raise ValueError("ジョブが完了していません。")
    result = None
    return_value = getattr(job, "return_value", None)
    if callable(return_value):
        result = return_value()
    else:
        result = return_value
    if result is None:
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore",
                category=DeprecationWarning,
                message="job.result is deprecated*",
            )
            result = job.result
    if result is None:
        raise ValueError("結果がまだ保存されていません。")
    return ZPEJobResultResponse(job_id=job.id, result=result)


@app.get("/calc/zpe/jobs/{job_id}/files")
def zpe_files(
    job_id: str,
    kind: str = Query("summary", pattern="^(summary|freqs|result)$"),
) -> FileResponse:
    try:
        job = fetch_job(job_id)
    except Exception as exc:  # pragma: no cover - depends on redis behavior
        raise HTTPException(status_code=404, detail="ジョブが見つかりません。") from exc
    status = job.get_status()
    if status == "failed":
        raise ValueError("ジョブが失敗しています。")
    if status != "finished":
        raise ValueError("ジョブが完了していません。")
    path = resolve_job_file(job, kind)
    return FileResponse(path, filename=path.name)
