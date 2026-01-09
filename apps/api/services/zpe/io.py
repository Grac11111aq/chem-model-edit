from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, Iterable

import numpy as np


def write_summary(
    path: Path,
    result: Dict[str, object],
    *,
    pseudo_dir: str,
    qe_input: str,
    new_calc: bool,
) -> None:
    lines = [
        "# ZPE summary (ASE Vibrations)",
        f"calc_start_time: {result['calc_start_time']}",
        f"calc_end_time: {result['calc_end_time']}",
        f"elapsed_seconds: {result['elapsed_seconds']:.1f}",
        f"qe_input: {qe_input}",
        f"pseudo_dir: {pseudo_dir}",
        f"selected_indices (mobile): {result['mobile_indices']}",
        f"delta (Å): {result['delta']}",
        f"kpts: {result['kpts']}",
        f"ecutwfc/ecutrho: {result.get('ecutwfc', 'N/A')}/{result.get('ecutrho', 'N/A')}",
        f"low_cut (cm^-1): {result['low_cut_cm']}",
        f"ZPE (eV): {result['zpe_ev']:.6f}",
        f"S_vib({result['temperature']} K): {result['s_vib_jmol_k']:.2f} J/mol/K",
        f"use_environ: {result['use_environ']}",
        f"calculation_mode: {'new' if new_calc else 'continue'}",
        f"cache_checked: {result['cache_checked']}",
        f"cache_deleted: {result['cache_deleted']}",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_freqs_csv(path: Path, freqs_cm: Iterable[float]) -> None:
    lines = ["frequency_cm^-1"]
    for freq in freqs_cm:
        lines.append(f"{float(np.real(freq)):.6f}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_result_json(path: Path, result: Dict[str, object]) -> None:
    path.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
