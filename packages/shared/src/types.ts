export type Atom = {
  symbol: string
  x: number
  y: number
  z: number
}

export type Structure = {
  atoms: Atom[]
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

export type Vector3 = {
  x: number
  y: number
  z: number
}

export type SupercellRequest = {
  structureA: Structure
  structureB: Structure
  sequence: string
  lattice: {
    a: Vector3
    b: Vector3
    c: Vector3
  }
}

export type SupercellResponse = {
  structure: Structure
  meta: {
    na: number
    nb: number
    layers: number
  }
}
