from __future__ import annotations

import math

from models import Lattice, LatticeParams, Vector3


def _length(vec: Vector3) -> float:
    return math.sqrt(vec.x**2 + vec.y**2 + vec.z**2)


def _dot(a: Vector3, b: Vector3) -> float:
    return a.x * b.x + a.y * b.y + a.z * b.z


def _clamp(value: float) -> float:
    return max(-1.0, min(1.0, value))


def vectors_to_params(lattice: Lattice) -> LatticeParams:
    a_len = _length(lattice.a)
    b_len = _length(lattice.b)
    c_len = _length(lattice.c)
    if a_len <= 0 or b_len <= 0 or c_len <= 0:
        raise ValueError("格子ベクトルの長さが無効です。")

    alpha = math.degrees(
        math.acos(_clamp(_dot(lattice.b, lattice.c) / (b_len * c_len)))
    )
    beta = math.degrees(math.acos(_clamp(_dot(lattice.a, lattice.c) / (a_len * c_len))))
    gamma = math.degrees(
        math.acos(_clamp(_dot(lattice.a, lattice.b) / (a_len * b_len)))
    )
    return LatticeParams(
        a=a_len,
        b=b_len,
        c=c_len,
        alpha=alpha,
        beta=beta,
        gamma=gamma,
    )


def params_to_vectors(params: LatticeParams) -> Lattice:
    a = params.a
    b = params.b
    c = params.c
    alpha = math.radians(params.alpha)
    beta = math.radians(params.beta)
    gamma = math.radians(params.gamma)

    if a <= 0 or b <= 0 or c <= 0:
        raise ValueError("格子定数が無効です。")

    sin_gamma = math.sin(gamma)
    if abs(sin_gamma) < 1.0e-8:
        raise ValueError("gamma が特異です。")

    a_vec = Vector3(x=a, y=0.0, z=0.0)
    b_vec = Vector3(x=b * math.cos(gamma), y=b * sin_gamma, z=0.0)

    c_x = c * math.cos(beta)
    c_y = c * (math.cos(alpha) - math.cos(beta) * math.cos(gamma)) / sin_gamma
    c_z_sq = c * c - c_x * c_x - c_y * c_y
    if c_z_sq < -1.0e-8:
        raise ValueError("格子角が不正で c_z を計算できません。")
    c_z = math.sqrt(max(c_z_sq, 0.0))
    c_vec = Vector3(x=c_x, y=c_y, z=c_z)

    return Lattice(a=a_vec, b=b_vec, c=c_vec)
