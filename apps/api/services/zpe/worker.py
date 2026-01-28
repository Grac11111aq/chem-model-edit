from __future__ import annotations

from contextlib import contextmanager
from dataclasses import dataclass
from datetime import datetime, timezone
import os
import shutil
import uuid
from pathlib import Path
from typing import Any, Dict

from ase.calculators.espresso import Espresso
from ase.vibrations import Vibrations
from rq import get_current_job

from models import ZPEJobRequest
from .cache import clean_vib_cache, sanitize_vib_cache
from .environ import prepare_environ_files, resolve_environ_path
from .io import format_freqs_csv, format_summary, write_result_json
from .parse import (
    ensure_mobile_indices,
    extract_fixed_indices,
    parse_atomic_species,
    parse_kpoints_automatic,
    parse_namelist,
    parse_qe_atoms,
)
from .paths import resolve_pseudo_dir, resolve_work_dir
from .qe import build_espresso_profile
from .result_store import get_result_store
from .settings import get_zpe_settings
from .thermo import calc_zpe_and_s_vib, normalize_frequencies


@contextmanager
def _chdir(path: Path):
    prev = Path.cwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(prev)


def _default_kpts() -> tuple[int, int, int]:
    return (1, 1, 1)


@dataclass
class ZPEComputeArtifacts:
    result: Dict[str, Any]
    summary_text: str
    freqs_csv: str


def _finalize_artifacts(
    *,
    job_dir: Path,
    result: Dict[str, Any],
    freqs_cm: list[float],
    pseudo_dir: str,
    qe_input_label: str,
    new_calc: bool,
) -> ZPEComputeArtifacts:
    summary_path = job_dir / "summary.txt"
    freqs_path = job_dir / "freqs.csv"
    result_path = job_dir / "result.json"

    summary_text = format_summary(
        result,
        pseudo_dir=pseudo_dir,
        qe_input=qe_input_label,
        new_calc=new_calc,
    )
    freqs_csv = format_freqs_csv(freqs_cm)

    summary_path.write_text(summary_text, encoding="utf-8")
    freqs_path.write_text(freqs_csv, encoding="utf-8")
    write_result_json(result_path, result)
    return ZPEComputeArtifacts(
        result=result,
        summary_text=summary_text,
        freqs_csv=freqs_csv,
    )


def _compute_mock_artifacts(
    request: ZPEJobRequest,
    *,
    settings: Any,
    job_dir: Path,
) -> ZPEComputeArtifacts:
    fixed_indices = extract_fixed_indices(request.content)
    kpts = parse_kpoints_automatic(request.content)
    kpts_use = kpts[0] if kpts else _default_kpts()

    mobile_indices = ensure_mobile_indices(
        request.mobile_indices,
        natoms=len(parse_qe_atoms(request.content)),
        fixed_indices=fixed_indices,
    )
    nfreq = max(1, len(mobile_indices) * 3)
    freqs_cm = [100.0 + 5.0 * i for i in range(nfreq)]
    zpe_ev, s_vib_jmol_k = calc_zpe_and_s_vib(
        freqs_cm,
        low_cut_cm=settings.low_cut_cm,
        temperature=settings.temperature,
    )

    now = datetime.now(timezone.utc).isoformat()
    result: Dict[str, Any] = {
        "freqs_cm": normalize_frequencies(freqs_cm),
        "zpe_ev": zpe_ev,
        "s_vib_jmol_k": s_vib_jmol_k,
        "mobile_indices": mobile_indices,
        "fixed_indices": fixed_indices,
        "kpts": kpts_use,
        "delta": settings.delta,
        "low_cut_cm": settings.low_cut_cm,
        "temperature": settings.temperature,
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
    return _finalize_artifacts(
        job_dir=job_dir,
        result=result,
        freqs_cm=freqs_cm,
        pseudo_dir="mock",
        qe_input_label="mock",
        new_calc=request.calc_mode == "new",
    )


def compute_zpe_artifacts(payload: Dict[str, Any], *, job_id: str) -> ZPEComputeArtifacts:
    settings = get_zpe_settings()
    request = ZPEJobRequest(**payload)

    work_dir = resolve_work_dir(settings)
    job_id_path = Path(job_id)
    if job_id_path.is_absolute() or job_id_path.parent != Path(".") or job_id_path.name in {"", ".", ".."}:
        raise ValueError(f"invalid job_id: {job_id}")
    job_dir = work_dir / job_id_path
    job_dir.mkdir(parents=True, exist_ok=True)

    if settings.worker_mode == "mock":
        return _compute_mock_artifacts(request, settings=settings, job_dir=job_dir)

    start_time = datetime.now(timezone.utc)
    with _chdir(job_dir):
        atoms = parse_qe_atoms(request.content)
        fixed_indices = extract_fixed_indices(request.content)
        mobile_indices = ensure_mobile_indices(
            request.mobile_indices,
            natoms=len(atoms),
            fixed_indices=fixed_indices,
        )

        pseudos = parse_atomic_species(request.content)
        if not pseudos:
            raise ValueError("ATOMIC_SPECIES section not found in QE input")

        kpts = parse_kpoints_automatic(request.content)
        kpts_use = kpts[0] if kpts else _default_kpts()

        control_defaults = {
            "calculation": "scf",
            "prefix": "ase_calc",
            "outdir": settings.outdir_name,
            "tprnfor": True,
            "tstress": False,
        }
        system_defaults = {
            "occupations": "smearing",
            "smearing": "mp",
            "degauss": 0.02,
        }
        electrons_defaults = {
            "conv_thr": 1.0e-6,
            "mixing_beta": 0.3,
            "mixing_mode": "plain",
            "electron_maxstep": 100,
        }

        control_from_in = parse_namelist(request.content, "control")
        system_from_in = parse_namelist(request.content, "system")
        electrons_from_in = parse_namelist(request.content, "electrons")

        pseudo_dir = resolve_pseudo_dir(
            request.content,
            request.input_dir,
            settings=settings,
        )

        control_dict = {**control_defaults, **control_from_in}
        control_dict["pseudo_dir"] = pseudo_dir

        system_dict = {**system_defaults, **system_from_in}
        electrons_dict = {**electrons_defaults, **electrons_from_in}

        input_data = {
            "control": control_dict,
            "system": system_dict,
            "electrons": electrons_dict,
        }

        calc_dir = job_dir / settings.calc_dir_name
        calc_dir.mkdir(parents=True, exist_ok=True)

        if request.use_environ:
            environ_src = resolve_environ_path(settings, request.input_dir, job_dir)
            prepare_environ_files(
                settings.vib_name,
                len(mobile_indices),
                environ_src,
                base_dir=job_dir,
            )
            shutil.copy2(environ_src, calc_dir / "environ.in")

        profile = build_espresso_profile(
            pseudo_dir=pseudo_dir,
            use_mpi=settings.use_mpi,
            np_core=settings.np_core,
            environ=request.use_environ,
            settings=settings,
        )

        atoms.calc = Espresso(
            pseudopotentials=pseudos,
            input_data=input_data,
            kpts=kpts_use,
            profile=profile,
            directory=str(calc_dir),
        )

        if request.calc_mode == "new":
            clean_vib_cache(job_dir, settings.vib_name)

        cache_info = sanitize_vib_cache(job_dir, settings.vib_name, natoms=len(atoms))

        vib = Vibrations(
            atoms,
            indices=mobile_indices,
            delta=settings.delta,
            name=settings.vib_name,
        )
        vib.run()
        freqs_cm = vib.get_frequencies()

    zpe_ev, s_vib_jmol_k = calc_zpe_and_s_vib(
        freqs_cm,
        low_cut_cm=settings.low_cut_cm,
        temperature=settings.temperature,
    )

    end_time = datetime.now(timezone.utc)
    elapsed_seconds = (end_time - start_time).total_seconds()

    qe_input_label = request.input_dir or "inline"
    checked_raw = cache_info.get("checked_files", 0)
    deleted_raw = cache_info.get("deleted_files", 0)
    checked_files = int(checked_raw) if isinstance(checked_raw, (int, float)) else 0
    deleted_files = int(deleted_raw) if isinstance(deleted_raw, (int, float)) else 0
    result: Dict[str, Any] = {
        "freqs_cm": normalize_frequencies(freqs_cm),
        "zpe_ev": zpe_ev,
        "s_vib_jmol_k": s_vib_jmol_k,
        "mobile_indices": mobile_indices,
        "fixed_indices": fixed_indices,
        "kpts": kpts_use,
        "delta": settings.delta,
        "low_cut_cm": settings.low_cut_cm,
        "temperature": settings.temperature,
        "use_environ": request.use_environ,
        "qe_input": qe_input_label,
        "pseudo_dir": pseudo_dir,
        "calc_start_time": start_time.isoformat(),
        "calc_end_time": end_time.isoformat(),
        "elapsed_seconds": elapsed_seconds,
        "cache_checked": checked_files,
        "cache_deleted": deleted_files,
        "ecutwfc": system_dict.get("ecutwfc"),
        "ecutrho": system_dict.get("ecutrho"),
    }

    return _finalize_artifacts(
        job_dir=job_dir,
        result=result,
        freqs_cm=freqs_cm,
        pseudo_dir=pseudo_dir,
        qe_input_label=qe_input_label,
        new_calc=request.calc_mode == "new",
    )


def run_zpe_job(payload: Dict[str, Any]) -> Dict[str, Any]:
    settings = get_zpe_settings()
    store = get_result_store()

    job = get_current_job()
    job_id = job.id if job else payload.get("job_id", f"local-{uuid.uuid4().hex[:8]}")
    store.set_status(job_id, "started")

    try:
        artifacts = compute_zpe_artifacts(payload, job_id=job_id)
        store.set_result(
            job_id,
            artifacts.result,
            summary_text=artifacts.summary_text,
            freqs_csv=artifacts.freqs_csv,
        )
        store.set_status(job_id, "finished")
        return artifacts.result
    except Exception as exc:
        store.set_status(job_id, "failed", detail=str(exc))
        raise
