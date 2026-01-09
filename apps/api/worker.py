from __future__ import annotations

from rq import Worker

from services.zpe import get_queue, get_redis_connection


def main() -> None:
    connection = get_redis_connection()
    queue = get_queue()
    worker = Worker([queue], connection=connection)
    worker.work(with_scheduler=True)


if __name__ == "__main__":
    main()
