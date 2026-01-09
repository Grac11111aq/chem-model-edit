from __future__ import annotations

import shutil
from pathlib import Path
from typing import Optional

from .settings import ZPESettings


def resolve_environ_path(
    settings: ZPESettings,
    input_dir: Optional[str],
    job_dir: Path,
) -> Path:
    if settings.environ_path:
        path = Path(settings.environ_path)
        if path.exists():
            return path
    candidates = []
    if input_dir:
        candidates.append(Path(input_dir) / "environ.in")
    candidates.append(job_dir / "environ.in")
    for candidate in candidates:
        if candidate.exists():
            return candidate
    raise FileNotFoundError("environ.in が見つかりません。ZPE_ENVIRON_PATH を指定してください。")


def prepare_environ_files(
    vib_name: str,
    n_mobile: int,
    environ_src: Path,
    *,
    base_dir: Path,
) -> None:
    eq_dir = base_dir / f"{vib_name}.eq"
    eq_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy2(environ_src, eq_dir / "environ.in")

    for i in range(n_mobile):
        for axis in ("x", "y", "z"):
            for sign in ("+", "-"):
                d = base_dir / f"{vib_name}.{i}{axis}{sign}"
                d.mkdir(parents=True, exist_ok=True)
                shutil.copy2(environ_src, d / "environ.in")

    alt_root = base_dir / "vib_ads_vac"
    if alt_root.exists():
        (alt_root / f"{vib_name}.eq").mkdir(parents=True, exist_ok=True)
        shutil.copy2(environ_src, alt_root / f"{vib_name}.eq" / "environ.in")
        for i in range(n_mobile):
            for axis in ("x", "y", "z"):
                for sign in ("+", "-"):
                    d = alt_root / f"{vib_name}.{i}{axis}{sign}"
                    d.mkdir(parents=True, exist_ok=True)
                    shutil.copy2(environ_src, d / "environ.in")
