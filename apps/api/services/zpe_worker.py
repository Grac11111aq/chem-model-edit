from __future__ import annotations

from contextlib import contextmanager
from datetime import datetime
import os
from pathlib import Path
import shutil
from typing import Any, Dict, Iterable, cast

from ase.constraints import FixAtoms
from ase.calculators.espresso import Espresso
from ase.vibrations import Vibrations
from rq import get_current_job

from services.zpe.cache import clean_vib_cache, sanitize_vib_cache
from services.zpe.environ import prepare_environ_files, resolve_environ_path
from services.zpe.io import write_freqs_csv, write_result_json, write_summary
from services.zpe.paths import resolve_pseudo_dir, resolve_work_dir
from services.zpe.parse import (
    ensure_mobile_indices,
    extract_fixed_indices,
    parse_atomic_species,
    parse_kpoints_automatic,
    parse_namelist,
    parse_qe_atoms,
)
from services.zpe.qe import build_espresso_profile
from services.zpe.settings import get_zpe_settings
from services.zpe.thermo import calc_zpe_and_s_vib, normalize_frequencies


@contextmanager
def _pushd(path: Path):
    prev = Path.cwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(prev)


def run_zpe_job(payload: Dict[str, Any]) -> Dict[str, object]:
    settings = get_zpe_settings()
    content = str(payload.get("content", ""))
    if not content.strip():
        raise ValueError("QE入力が空です。")

    work_dir = resolve_work_dir(settings)
    job = get_current_job()
    job_id = job.id if job else datetime.now().strftime("%Y%m%d%H%M%S")
    job_dir = work_dir / job_id
    job_dir.mkdir(parents=True, exist_ok=True)

    input_path = job_dir / "input.in"
    input_path.write_text(content, encoding="utf-8")

    fixed_indices = extract_fixed_indices(content)
    atoms = parse_qe_atoms(content)
    if fixed_indices:
        atoms.set_constraint(FixAtoms(indices=fixed_indices))

    raw_mobile = payload.get("mobile_indices", [])
    mobile_indices = ensure_mobile_indices(
        cast(Iterable[int], raw_mobile),
        len(atoms),
        fixed_indices,
    )

    input_dir = payload.get("input_dir")
    input_dir_str = str(input_dir) if isinstance(input_dir, str) and input_dir else None
    pseudo_dir = resolve_pseudo_dir(content, input_dir_str, settings=settings)

    use_environ = bool(payload.get("use_environ", False))
    new_calc = bool(payload.get("new_calc", False))

    temperature_raw = payload.get("temperature")
    low_cut_raw = payload.get("low_cut_cm")
    temperature = float(temperature_raw) if temperature_raw is not None else settings.temperature
    low_cut_cm = float(low_cut_raw) if low_cut_raw is not None else settings.low_cut_cm
    delta = settings.delta

    system_from_in = parse_namelist(content, "system")
    electrons_from_in = parse_namelist(content, "electrons")

    sys_defaults = {
        "occupations": "smearing",
        "smearing": "mp",
        "degauss": 0.02,
    }
    ele_defaults = {
        "conv_thr": 1.0e-6,
        "mixing_beta": 0.3,
        "mixing_mode": "plain",
        "electron_maxstep": 100,
    }
    system_dict = {**sys_defaults, **system_from_in}
    electrons_dict = {**ele_defaults, **electrons_from_in}

    kpoints = parse_kpoints_automatic(content)
    kpts_use = (3, 3, 1)
    if kpoints is not None:
        (nkx, nky, nkz), _ = kpoints
        kpts_use = (nkx, nky, nkz)

    pseudos = parse_atomic_species(content)
    if not pseudos:
        raise ValueError("ATOMIC_SPECIES が見つからず pseudopotentials を構成できません。")

    calc_dir = job_dir / settings.calc_dir_name
    calc_dir.mkdir(parents=True, exist_ok=True)

    if use_environ:
        environ_src = resolve_environ_path(settings, input_dir_str, job_dir)
        prepare_environ_files(
            settings.vib_name,
            len(mobile_indices),
            environ_src,
            base_dir=job_dir,
        )
        shutil.copy2(environ_src, calc_dir / "environ.in")

    input_data = {
        "control": {
            "calculation": "scf",
            "prefix": "ase_calc",
            "pseudo_dir": pseudo_dir,
            "outdir": str(job_dir / settings.outdir_name),
            "tprnfor": True,
            "tstress": False,
        },
        "system": system_dict,
        "electrons": electrons_dict,
    }

    profile = build_espresso_profile(
        pseudo_dir=pseudo_dir,
        use_mpi=settings.use_mpi,
        np_core=settings.np_core,
        environ=use_environ,
        settings=settings,
    )

    if new_calc:
        clean_vib_cache(job_dir, settings.vib_name)

    start_time = datetime.now()
    with _pushd(job_dir):
        cleaned = sanitize_vib_cache(job_dir, settings.vib_name, natoms=len(atoms))

        atoms.calc = Espresso(
            pseudopotentials=pseudos,
            input_data=input_data,
            kpts=kpts_use,
            profile=profile,
            directory=str(calc_dir),
        )
        vib = Vibrations(atoms, indices=mobile_indices, delta=delta, name=settings.vib_name)
        vib.run()
        freqs_cm = vib.get_frequencies()

    freqs_real = normalize_frequencies(freqs_cm)
    zpe_ev, s_vib_jmol_k = calc_zpe_and_s_vib(
        freqs_cm,
        low_cut_cm=low_cut_cm,
        temperature=temperature,
    )

    end_time = datetime.now()
    elapsed_seconds = (end_time - start_time).total_seconds()

    result = {
        "frequencies_cm": freqs_real,
        "zpe_ev": float(zpe_ev),
        "s_vib_jmol_k": float(s_vib_jmol_k),
        "low_cut_cm": float(low_cut_cm),
        "temperature": float(temperature),
        "delta": float(delta),
        "kpts": [int(kpts_use[0]), int(kpts_use[1]), int(kpts_use[2])],
        "ecutwfc": system_dict.get("ecutwfc"),
        "ecutrho": system_dict.get("ecutrho"),
        "mobile_indices": list(mobile_indices),
        "fixed_indices": list(fixed_indices),
        "use_environ": bool(use_environ),
        "calculation_mode": "new" if new_calc else "continue",
        "calc_start_time": start_time.isoformat(),
        "calc_end_time": end_time.isoformat(),
        "elapsed_seconds": float(elapsed_seconds),
        "cache_checked": int(cast(int, cleaned["checked_files"])),
        "cache_deleted": int(cast(int, cleaned["deleted_files"])),
    }

    summary_path = job_dir / "zpe_summary.txt"
    csv_path = job_dir / "zpe_freqs.csv"
    result_path = job_dir / "zpe_result.json"

    write_summary(
        summary_path,
        result,
        pseudo_dir=pseudo_dir,
        qe_input=str(input_path),
        new_calc=new_calc,
    )
    write_freqs_csv(csv_path, freqs_cm)
    write_result_json(result_path, result)

    if job is not None:
        job.meta.update(
            {
                "work_dir": str(job_dir.resolve()),
                "summary_file": str(summary_path.resolve()),
                "freqs_file": str(csv_path.resolve()),
                "result_file": str(result_path.resolve()),
            }
        )
        job.save_meta()

    return result
