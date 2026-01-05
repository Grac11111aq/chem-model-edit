export type Atom = {
  symbol: string
  x: number
  y: number
  z: number
}

export type Structure = {
  atoms: Atom[]
  lattice?: Lattice | null
}

export type Vector3 = {
  x: number
  y: number
  z: number
}

export type Lattice = {
  a: Vector3
  b: Vector3
  c: Vector3
}

export type SupercellMeta = {
  na: number
  nb: number
  layers: number
}
