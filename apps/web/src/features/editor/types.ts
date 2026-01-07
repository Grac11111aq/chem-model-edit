import type { Atom, Lattice } from '@/lib/types'

export type AtomDraft = {
  symbol: string
  x: string
  y: string
  z: string
}

export type StructureState = {
  id: string
  name: string
  color: string
  atoms: Array<Atom>
  drafts: Array<AtomDraft>
  isVisible: boolean
  opacity: number
  lattice?: Lattice | null
}

export type Matrix3 = [
  [number, number, number],
  [number, number, number],
  [number, number, number],
]

export type PbcState =
  | { enabled: false; reason: string; mismatch?: boolean }
  | { enabled: true; reason: string; matrix: Matrix3; inverse: Matrix3 }
