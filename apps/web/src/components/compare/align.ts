import type { Atom } from '../../lib/types'

type Vec3 = { x: number; y: number; z: number }

export function shiftAtoms(
  atoms: Array<Atom>,
  indices: Array<number>,
  delta: Vec3,
): Array<Atom> {
  if (indices.length === 0) return atoms
  const set = new Set(indices)
  return atoms.map((atom, index) =>
    set.has(index)
      ? {
          ...atom,
          x: atom.x + delta.x,
          y: atom.y + delta.y,
          z: atom.z + delta.z,
        }
      : atom,
  )
}

export function alignSelectedToOrigin(
  atoms: Array<Atom>,
  indices: Array<number>,
): Array<Atom> {
  if (indices.length === 0) return atoms
  const anchorIndex = indices[0]
  if (anchorIndex < 0 || anchorIndex >= atoms.length) return atoms
  const anchor = atoms[anchorIndex]
  return shiftAtoms(atoms, indices, {
    x: -anchor.x,
    y: -anchor.y,
    z: -anchor.z,
  })
}

export function alignSelectedCentroid(
  atoms: Array<Atom>,
  indices: Array<number>,
): Array<Atom> {
  if (indices.length === 0) return atoms
  const selectedIndices = indices.filter(
    (index) => index >= 0 && index < atoms.length,
  )
  if (selectedIndices.length === 0) return atoms
  const selected = selectedIndices.map((index) => atoms[index])
  const sum = selected.reduce(
    (acc, atom) => ({
      x: acc.x + atom.x,
      y: acc.y + atom.y,
      z: acc.z + atom.z,
    }),
    { x: 0, y: 0, z: 0 },
  )
  const inv = 1 / selected.length
  return shiftAtoms(atoms, indices, {
    x: -sum.x * inv,
    y: -sum.y * inv,
    z: -sum.z * inv,
  })
}
