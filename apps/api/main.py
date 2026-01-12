from __future__ import annotations

import logging
import secrets
from typing import Any, Literal

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
    StructureCreateRequest,
    StructureCreateResponse,
    SupercellBuildMeta,
    SupercellBuildRequest,
    SupercellBuildResponse,
    SupercellRequest,
    SupercellResponse,
    TiledSupercellRequest,
    ZPEJobRequest,
    ZPEJobResponse,
    ZPEJobResultResponse,
    ZPEJobStatusResponse,
    ZPEComputeRegisterRequest,
    ZPEComputeRegisterResponse,
    ZPEEnrollTokenRequest,
    ZPEEnrollTokenResponse,
    ZPEParseRequest,
    ZPEParseResponse,
    ZPEResult,
)
from services.export import export_qe_in
from services.lattice import params_to_vectors, vectors_to_params
from services.parse import parse_qe_in, structure_from_ase
from services.structures import (
    create_structure_from_qe,
    get_structure_bcif,
    get_structure_entry,
    register_structure_atoms,
)
from services.supercell import (
    build_supercell_from_grid,
    generate_supercell,
    generate_tiled_supercell,
)
from services.transplant import transplant_delta
from services.zpe import (
    ensure_mobile_indices,
    enqueue_zpe_job,
    get_result_store,
    get_zpe_settings,
)
from services.zpe.enroll import get_enroll_store
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


def _cell_has_lattice(cell: Any) -> bool:
    try:
        lengths = cell.lengths()
    except AttributeError:
        return False
    return any(length > 1.0e-8 for length in lengths)


def _cells_match(cell_a: Any, cell_b: Any, *, tol: float = 1.0e-6) -> bool:
    try:
        arr_a = cell_a.array
        arr_b = cell_b.array
    except AttributeError:
        return False
    for row_a, row_b in zip(arr_a, arr_b, strict=True):
        for value_a, value_b in zip(row_a, row_b, strict=True):
            if abs(float(value_a) - float(value_b)) > tol:
                return False
    return True


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


def _extract_admin_token(request: Request) -> str | None:
    token: str | None
    auth = request.headers.get("authorization")
    if auth and auth.lower().startswith("bearer "):
        token = auth.split(" ", 1)[1].strip()
        return token or None
    token = request.headers.get("x-admin-token")
    return token or None


def require_admin(request: Request) -> None:
    settings = get_zpe_settings()
    if not settings.admin_token:
        return
    token = _extract_admin_token(request)
    if not token or not secrets.compare_digest(token, settings.admin_token):
        raise HTTPException(status_code=401, detail="unauthorized")


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


@app.post("/structures", response_model=StructureCreateResponse)
def create_structure(request: StructureCreateRequest) -> StructureCreateResponse:
    structure_id, structure, source = create_structure_from_qe(request.content)
    return StructureCreateResponse(
        structure_id=structure_id,
        structure=structure,
        source=source,
    )


@app.get("/structures/{structure_id}/view")
def view_structure(
    structure_id: str,
    format: str = Query("bcif"),
    lossy: int = Query(0, ge=0, le=1),
    precision: int = Query(3, ge=0),
) -> Response:
    if format != "bcif":
        raise HTTPException(status_code=400, detail="Unsupported format")

    try:
        bcif = get_structure_bcif(
            structure_id,
            lossy=lossy == 1,
            precision=precision if lossy == 1 else None,
        )
    except KeyError as exc:
        raise HTTPException(status_code=404, detail="Structure not found") from exc

    return Response(content=bcif, media_type="application/x-bcif")


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


@app.post("/supercell/build", response_model=SupercellBuildResponse)
def supercell_build(request: SupercellBuildRequest) -> SupercellBuildResponse:
    try:
        base_entry = get_structure_entry(request.baseStructureId)
    except KeyError as exc:
        raise HTTPException(status_code=404, detail="Structure not found") from exc

    base_atoms = base_entry.atoms
    base_cell = base_atoms.get_cell()
    if not _cell_has_lattice(base_cell):
        raise HTTPException(status_code=400, detail="base structure has no lattice")

    tiles: dict[str, Any] = {}
    structure_ids_used: list[str] = []
    for row in request.grid.tiles:
        for structure_id in row:
            if structure_id in tiles:
                continue
            structure_ids_used.append(structure_id)
            try:
                entry = get_structure_entry(structure_id)
            except KeyError as exc:
                raise HTTPException(
                    status_code=404, detail="Structure not found"
                ) from exc
            tiles[structure_id] = entry.atoms

    options = request.options
    if options and options.validateLattice != "none":
        for structure_id, atoms in tiles.items():
            if structure_id == request.baseStructureId:
                continue
            cell = atoms.get_cell()
            lattice_ok = _cell_has_lattice(cell) and _cells_match(base_cell, cell)
            if lattice_ok:
                continue
            if options.validateLattice == "error":
                raise HTTPException(
                    status_code=400,
                    detail=f"lattice mismatch for structure {structure_id}",
                )
            logger.warning(
                "supercell.build lattice mismatch: %s vs base %s",
                structure_id,
                request.baseStructureId,
            )

    check_overlap = bool(options.checkOverlap) if options else False
    overlap_tolerance = options.overlapTolerance if options else None
    atoms_out, overlap_count = build_supercell_from_grid(
        base_atoms,
        request.grid,
        tiles,
        check_overlap=check_overlap,
        overlap_tolerance=overlap_tolerance,
    )

    structure_id = register_structure_atoms(atoms_out, source="supercell-build")
    include_structure = (
        bool(request.output.includeStructure) if request.output else False
    )
    structure = structure_from_ase(atoms_out) if include_structure else None
    meta = SupercellBuildMeta(
        rows=request.grid.rows,
        cols=request.grid.cols,
        tileCount=request.grid.rows * request.grid.cols,
        overlapCount=overlap_count if check_overlap else None,
        baseStructureId=request.baseStructureId,
        structureIdsUsed=structure_ids_used,
    )
    return SupercellBuildResponse(
        structureId=structure_id,
        structure=structure,
        meta=meta,
    )


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

@app.post("/calc/zpe/compute/enroll-tokens", response_model=ZPEEnrollTokenResponse)
def zpe_compute_enroll_token(
    request: ZPEEnrollTokenRequest,
    raw: Request,
) -> ZPEEnrollTokenResponse:
    require_admin(raw)
    if request.ttl_seconds is not None and request.ttl_seconds <= 0:
        raise HTTPException(status_code=400, detail="ttl_seconds must be >= 1")
    store = get_enroll_store()
    token = store.create_token(ttl_seconds=request.ttl_seconds, label=request.label)
    return ZPEEnrollTokenResponse(
        token=token.token,
        expires_at=token.expires_at,
        ttl_seconds=token.ttl_seconds,
        label=token.label,
    )


@app.post("/calc/zpe/compute/servers/register", response_model=ZPEComputeRegisterResponse)
def zpe_compute_register(
    request: ZPEComputeRegisterRequest,
) -> ZPEComputeRegisterResponse:
    store = get_enroll_store()
    try:
        registration = store.consume_token(
            request.token,
            name=request.name,
            meta=request.meta,
        )
    except KeyError as exc:
        raise HTTPException(status_code=400, detail="invalid enroll token") from exc
    return ZPEComputeRegisterResponse(
        server_id=registration.server_id,
        registered_at=registration.registered_at,
        name=registration.name,
    )
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
