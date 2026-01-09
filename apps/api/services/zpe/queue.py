from __future__ import annotations

from redis import Redis
from rq import Queue
from rq.job import Job

from .settings import get_zpe_settings


def get_redis_connection() -> Redis:
    settings = get_zpe_settings()
    return Redis.from_url(settings.redis_url)


def get_queue() -> Queue:
    settings = get_zpe_settings()
    return Queue(settings.queue_name, connection=get_redis_connection())


def fetch_job(job_id: str) -> Job:
    return Job.fetch(job_id, connection=get_redis_connection())
