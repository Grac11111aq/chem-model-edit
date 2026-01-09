from __future__ import annotations

from functools import lru_cache
from typing import Optional

from pydantic_settings import BaseSettings, SettingsConfigDict


class ZPESettings(BaseSettings):
    model_config = SettingsConfigDict(
        env_prefix="ZPE_",
        env_file=".env",
        extra="ignore",
    )

    redis_url: str = "redis://localhost:6379/0"
    work_dir: str = "zpe_jobs"
    queue_name: str = "zpe"

    pw_command: str = "pw.x"
    pw_path: Optional[str] = None
    use_mpi: bool = True
    mpi_cmd: str = "mpirun"
    np_core: int = 12

    pseudo_dir: Optional[str] = None
    allow_input_pseudo_dir: bool = False

    delta: float = 0.01
    low_cut_cm: float = 50.0
    temperature: float = 298.15
    vib_name: str = "vib_ads"
    outdir_name: str = "calc_scratch"
    calc_dir_name: str = "vib_ads_vac"
    job_timeout_seconds: int = 86400
    result_ttl_seconds: int = 604800

    environ_path: Optional[str] = None


@lru_cache
def get_zpe_settings() -> ZPESettings:
    return ZPESettings()
