from __future__ import annotations

import json
import shutil
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np


def sanitize_vib_cache(
    job_dir: Path, name: str, natoms: Optional[int] = None
) -> Dict[str, object]:
    roots_relaxed = [job_dir / name]

    deleted = 0
    checked = 0
    deleted_list: List[str] = []

    def check_and_maybe_delete(jpath: Path, strict: bool) -> None:
        nonlocal deleted, checked, deleted_list
        checked += 1
        try:
            if not jpath.exists():
                return
            if jpath.stat().st_size == 0:
                raise ValueError("empty cache")
            data = json.loads(jpath.read_text(encoding="utf-8"))
            forces = data.get("forces")
            if forces is None or (isinstance(forces, list) and len(forces) == 0):
                raise ValueError("forces missing/empty")
            if strict:
                if natoms is not None and len(forces) != natoms:
                    raise ValueError("forces length mismatch")
                arr = np.array(forces, dtype=float)
                if arr.ndim != 2 or arr.shape[1] != 3:
                    raise ValueError("forces bad shape")
                if not np.isfinite(arr).all():
                    raise ValueError("forces contains NaN/inf")
        except Exception:
            try:
                jpath.unlink()
                deleted += 1
                deleted_list.append(str(jpath))
            except Exception:
                pass

    for root in (job_dir, job_dir / "vib_ads_vac"):
        eq_json = root / f"{name}.eq" / "cache.eq.json"
        if eq_json.exists():
            check_and_maybe_delete(eq_json, strict=True)
        for d in root.glob(f"{name}.*"):
            if d.is_dir():
                for jf in d.glob("cache*.json"):
                    check_and_maybe_delete(jf, strict=True)

    for base in roots_relaxed:
        if base.exists() and base.is_dir():
            for jf in base.rglob("cache*.json"):
                check_and_maybe_delete(jf, strict=False)

    for base in (job_dir / "vib_ads_vac" / name,):
        if base.exists() and base.is_dir():
            for jf in base.rglob("cache*.json"):
                check_and_maybe_delete(jf, strict=True)

    return {
        "deleted_files": deleted,
        "checked_files": checked,
        "deleted_list": deleted_list,
    }


def clean_vib_cache(job_dir: Path, vib_name: str) -> None:
    for path in job_dir.glob(f"{vib_name}.*"):
        shutil.rmtree(path, ignore_errors=True)
