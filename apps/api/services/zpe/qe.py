from __future__ import annotations

import inspect
import shutil
from pathlib import Path
from typing import Any, Dict, List

from ase.calculators.espresso import EspressoProfile

from .settings import ZPESettings


def build_espresso_profile(
    *,
    pseudo_dir: str,
    use_mpi: bool,
    np_core: int,
    environ: bool,
    settings: ZPESettings,
) -> EspressoProfile:
    pw = resolve_pw_command(settings)
    argv: List[str] = []
    if use_mpi:
        argv += [settings.mpi_cmd, "-np", str(np_core), pw]
    else:
        argv += [pw]
    if environ:
        argv += ["-environ"]

    kwargs: Dict[str, Any] = {}
    sig = inspect.signature(EspressoProfile)
    params = set(sig.parameters)
    if "pseudo_dir" in params:
        kwargs["pseudo_dir"] = pseudo_dir

    if {"command", "argv"} <= params:
        kwargs["command"] = settings.mpi_cmd if use_mpi else pw
        if use_mpi:
            kwargs["argv"] = ["-np", str(np_core), pw] + (
                ["-environ"] if environ else []
            )
        else:
            kwargs["argv"] = ["-environ"] if environ else []
    elif "command" in params:
        kwargs["command"] = " ".join(argv)
    elif "argv" in params:
        kwargs["argv"] = argv
    else:
        raise RuntimeError("Unsupported EspressoProfile signature")
    return EspressoProfile(**kwargs)


def resolve_pw_command(settings: ZPESettings) -> str:
    if settings.pw_path:
        path = Path(settings.pw_path)
        if path.exists():
            return str(path)
    pw = shutil.which(settings.pw_command)
    if pw is None:
        raise FileNotFoundError(
            "pw.x が見つかりません。ZPE_PW_PATH を設定してください。"
        )
    return pw
