import type { Lattice, Vector3 } from '@/lib/types'
import type { Matrix3 } from './types'

export const LATTICE_TOLERANCE = 1.0e-3

export const latticeToMatrix = (lattice: Lattice): Matrix3 => [
  [lattice.a.x, lattice.a.y, lattice.a.z],
  [lattice.b.x, lattice.b.y, lattice.b.z],
  [lattice.c.x, lattice.c.y, lattice.c.z],
]

export const maxAbsDifference = (a: Matrix3, b: Matrix3) => {
  let max = 0
  for (let i = 0; i < 3; i += 1) {
    for (let j = 0; j < 3; j += 1) {
      const diff = Math.abs(a[i][j] - b[i][j])
      if (diff > max) {
        max = diff
      }
    }
  }
  return max
}

export const invert3 = (m: Matrix3): Matrix3 | null => {
  const [[a00, a01, a02], [a10, a11, a12], [a20, a21, a22]] = m
  const b01 = a22 * a11 - a12 * a21
  const b11 = -a22 * a10 + a12 * a20
  const b21 = a21 * a10 - a11 * a20
  const det = a00 * b01 + a01 * b11 + a02 * b21
  if (Math.abs(det) < 1.0e-12) {
    return null
  }
  const invDet = 1 / det
  return [
    [
      b01 * invDet,
      (-a22 * a01 + a02 * a21) * invDet,
      (a12 * a01 - a02 * a11) * invDet,
    ],
    [
      b11 * invDet,
      (a22 * a00 - a02 * a20) * invDet,
      (-a12 * a00 + a02 * a10) * invDet,
    ],
    [
      b21 * invDet,
      (-a21 * a00 + a01 * a20) * invDet,
      (a11 * a00 - a01 * a10) * invDet,
    ],
  ]
}

export const multiplyMatVec = (m: Matrix3, v: Vector3): Vector3 => ({
  x: m[0][0] * v.x + m[0][1] * v.y + m[0][2] * v.z,
  y: m[1][0] * v.x + m[1][1] * v.y + m[1][2] * v.z,
  z: m[2][0] * v.x + m[2][1] * v.y + m[2][2] * v.z,
})

export const wrapFractional = (v: Vector3): Vector3 => ({
  x: v.x - Math.round(v.x),
  y: v.y - Math.round(v.y),
  z: v.z - Math.round(v.z),
})

export const minimumImageDelta = (
  delta: Vector3,
  lattice: Matrix3,
  inverse: Matrix3,
): Vector3 => {
  const fractional = multiplyMatVec(inverse, delta)
  const wrapped = wrapFractional(fractional)
  return multiplyMatVec(lattice, wrapped)
}
