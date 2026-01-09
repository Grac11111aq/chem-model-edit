from __future__ import annotations

import math
from typing import Iterable, List, Tuple

import numpy as np
from ase.units import kB

J_PER_EV = 1.602176634e-19
H_PLANCK_J = 6.62607015e-34  # J*s (CODATA 2023)
H_PLANCK_EV = H_PLANCK_J / J_PER_EV
HC_EV_CM = 1.239841984e-4
C_LIGHT = 2.99792458e10


def normalize_frequencies(freqs_cm: Iterable[float]) -> List[float]:
    freqs_cm_arr = np.asarray(list(freqs_cm), dtype=np.complex128)
    return [float(np.real(v)) for v in freqs_cm_arr]


def _select_valid_frequencies(
    freqs_cm: Iterable[float],
    *,
    low_cut_cm: float,
    imag_tol: float = 1.0e-6,
) -> List[float]:
    freqs_cm_arr = np.asarray(list(freqs_cm), dtype=np.complex128)
    freqs_real = np.real(freqs_cm_arr)
    freqs_imag = np.imag(freqs_cm_arr)
    valid: List[float] = []
    for f_real, f_imag in zip(freqs_real, freqs_imag):
        if not np.isfinite(f_real):
            continue
        if abs(f_imag) > imag_tol:
            continue
        if f_real <= low_cut_cm:
            continue
        valid.append(float(f_real))
    return valid


def calc_zpe_and_s_vib(
    freqs_cm: Iterable[float],
    *,
    low_cut_cm: float,
    temperature: float,
) -> Tuple[float, float]:
    valid = _select_valid_frequencies(freqs_cm, low_cut_cm=low_cut_cm)
    zpe_ev = 0.5 * HC_EV_CM * float(np.sum(valid)) if valid else 0.0

    freqs_hz = np.array(valid, dtype=float) * C_LIGHT
    s_vib = 0.0
    for nu in freqs_hz:
        x = H_PLANCK_EV * nu / (kB * temperature)
        if not math.isfinite(x) or x <= 1.0e-12:
            continue
        exp_x = math.exp(x)
        denom = exp_x - 1.0
        if denom <= 0.0:
            continue
        exp_neg = math.exp(-x)
        if exp_neg >= 1.0:
            continue
        s_vib += kB * (x / denom - math.log1p(-exp_neg))
    s_vib_jmol_k = float(s_vib * 96485.33212)
    return float(zpe_ev), s_vib_jmol_k
