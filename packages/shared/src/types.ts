export type Atom = {
  symbol: string
  x: number
  y: number
  z: number
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

export type Structure = {
  atoms: Atom[]
  lattice?: Lattice | null
}

export type ParseRequest = {
  content: string
  format?: string
}

export type ParseResponse = {
  structure: Structure
}

export type ExportRequest = {
  structure: Structure
  format?: string
}

export type ExportResponse = {
  content: string
}

export type SupercellRequest = {
  structureA: Structure
  structureB: Structure
  sequence: string
  lattice: Lattice
}

export type SupercellResponse = {
  structure: Structure
  meta: {
    na: number
    nb: number
    layers: number
  }
}
