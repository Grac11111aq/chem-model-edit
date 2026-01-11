from __future__ import annotations

from typing import Any, Dict

from .queue import get_queue
from .result_store import get_result_store
from .settings import get_zpe_settings


def enqueue_zpe_job(payload: Dict[str, Any]) -> str:
    settings = get_zpe_settings()
    if settings.compute_mode != "remote-queue":
        raise ValueError("compute_mode must be 'remote-queue'.")
    store = get_result_store()
    job = get_queue().enqueue(
        "services.zpe.worker.run_zpe_job",
        payload,
        job_timeout=settings.job_timeout_seconds,
        result_ttl=settings.result_ttl_seconds,
    )
    store.set_status(job.id, "queued")
    return job.id
