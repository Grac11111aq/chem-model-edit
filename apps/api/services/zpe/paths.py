from __future__ import annotations

from pathlib import Path
from typing import Optional

from rq.job import Job

from .parse import parse_namelist
from .settings import ZPESettings, get_zpe_settings


def resolve_pseudo_dir(
    content: str,
    input_dir: Optional[str],
    *,
    settings: Optional[ZPESettings] = None,
) -> str:
    settings = settings or get_zpe_settings()
    if settings.pseudo_dir:
        return settings.pseudo_dir
    if settings.allow_input_pseudo_dir:
        control = parse_namelist(content, "control")
        pseudo_dir = control.get("pseudo_dir")
        if isinstance(pseudo_dir, str) and pseudo_dir:
            path = Path(pseudo_dir)
            if path.is_absolute():
                return str(path)
            if input_dir:
                return str(Path(input_dir) / path)
    raise ValueError("pseudo_dir が設定されていません。ZPE_PSEUDO_DIR を設定してください。")


def resolve_work_dir(settings: Optional[ZPESettings] = None) -> Path:
    settings = settings or get_zpe_settings()
    base = Path(settings.work_dir).expanduser()
    if not base.is_absolute():
        base = (Path.cwd() / base).resolve()
    else:
        base = base.resolve()
    base.mkdir(parents=True, exist_ok=True)
    return base


def resolve_job_file(job: Job, kind: str) -> Path:
    key_map = {
        "summary": "summary_file",
        "freqs": "freqs_file",
        "result": "result_file",
    }
    key = key_map.get(kind)
    if key is None:
        raise ValueError("kind は summary / freqs / result のいずれかです。")
    path_value = job.meta.get(key)
    if not path_value:
        raise ValueError("結果ファイルがまだ準備できていません。")
    path = Path(str(path_value)).resolve()
    work_dir_value = job.meta.get("work_dir")
    if work_dir_value:
        work_dir = Path(str(work_dir_value)).expanduser().resolve()
    else:
        work_dir = resolve_work_dir().resolve()
    if not path.exists():
        raise ValueError("結果ファイルが見つかりません。")
    if not path.is_relative_to(work_dir):
        raise ValueError("結果ファイルのパスが不正です。")
    return path
