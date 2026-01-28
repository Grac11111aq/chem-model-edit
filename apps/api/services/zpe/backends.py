from __future__ import annotations

from datetime import datetime, timezone
from typing import Dict, Any
from uuid import uuid4

from models import ZPEJobRequest
from .io import format_freqs_csv, format_summary
from .parse import extract_fixed_indices, parse_kpoints_automatic
from .http_queue import enqueue_http_job
from .queue import get_queue
from .result_store import ResultStore, get_result_store
from .settings import get_zpe_settings
from .thermo import calc_zpe_and_s_vib, normalize_frequencies


def _now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def _default_kpts() -> tuple[int, int, int]:
    return (1, 1, 1)


def enqueue_zpe_job(payload: Dict[str, Any], *, queue_name: str | None = None) -> str:
    settings = get_zpe_settings()
    if settings.compute_mode not in {"remote-queue", "remote-http", "mock"}:
        raise ValueError(
            "compute_mode must be 'remote-queue', 'remote-http', or 'mock'."
        )
    store = get_result_store()
    if settings.compute_mode == "remote-queue":
        job = get_queue(queue_name).enqueue(
            "services.zpe.worker.run_zpe_job",
            payload,
            job_timeout=settings.job_timeout_seconds,
            result_ttl=settings.result_ttl_seconds,
        )
        store.set_status(job.id, "queued")
        return job.id
    if settings.compute_mode == "remote-http":
        return enqueue_http_job(payload)
    return _run_mock_job(payload, store)


def _run_mock_job(payload: Dict[str, Any], store: ResultStore) -> str:
    request = ZPEJobRequest(**payload)
    job_id = f"mock-{uuid4().hex}"
    store.set_status(job_id, "started")

    fixed_indices = extract_fixed_indices(request.content)
    kpts = parse_kpoints_automatic(request.content)
    kpts_use = kpts[0] if kpts else _default_kpts()

    nfreq = max(1, len(request.mobile_indices) * 3)
    freqs_cm = [100.0 + 5.0 * i for i in range(nfreq)]
    zpe_ev, s_vib_jmol_k = calc_zpe_and_s_vib(
        freqs_cm,
        low_cut_cm=get_zpe_settings().low_cut_cm,
        temperature=get_zpe_settings().temperature,
    )

    now = _now_iso()
    result: Dict[str, Any] = {
        "freqs_cm": normalize_frequencies(freqs_cm),
        "zpe_ev": zpe_ev,
        "s_vib_jmol_k": s_vib_jmol_k,
        "mobile_indices": request.mobile_indices,
        "fixed_indices": fixed_indices,
        "kpts": kpts_use,
        "delta": get_zpe_settings().delta,
        "low_cut_cm": get_zpe_settings().low_cut_cm,
        "temperature": get_zpe_settings().temperature,
        "use_environ": request.use_environ,
        "qe_input": "mock",
        "pseudo_dir": "mock",
        "calc_start_time": now,
        "calc_end_time": now,
        "elapsed_seconds": 0.0,
        "cache_checked": 0,
        "cache_deleted": 0,
        "ecutwfc": None,
        "ecutrho": None,
    }

    summary_text = format_summary(
        result,
        pseudo_dir="mock",
        qe_input="mock",
        new_calc=request.calc_mode == "new",
    )
    freqs_csv = format_freqs_csv(freqs_cm)
    store.set_result(job_id, result, summary_text=summary_text, freqs_csv=freqs_csv)
    store.set_status(job_id, "finished")
    return job_id
