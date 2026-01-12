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

export type LatticeParams = {
  a: number
  b: number
  c: number
  alpha: number
  beta: number
  gamma: number
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

export type StructureRegisterRequest = {
  structure: Structure
  source?: string
}

export type StructureRegisterResponse = {
  structureId: string
  source: string
}

export type StructureImportRequest = {
  content: string
  format?: string
}

export type StructureImportResponse = {
  structureId: string
  structure: Structure
  source: string
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
    overlapCount?: number
  }
}

export type SupercellGridAxis =
  | { row: 'a'; col: 'b' }
  | { row: 'b'; col: 'a' }

export type SupercellGrid = {
  rows: number
  cols: number
  tiles: string[][]
  axis?: SupercellGridAxis
}

export type SupercellBuildOptions = {
  checkOverlap?: boolean
  overlapTolerance?: number
  validateLattice?: 'none' | 'warn' | 'error'
}

export type SupercellBuildOutput = {
  includeStructure?: boolean
}

export type SupercellBuildRequest = {
  baseStructureId: string
  grid: SupercellGrid
  options?: SupercellBuildOptions
  output?: SupercellBuildOutput
}

export type SupercellBuildMeta = {
  rows: number
  cols: number
  tileCount: number
  overlapCount?: number
  baseStructureId: string
  structureIdsUsed: string[]
}

export type SupercellBuildResponse = {
  structureId: string
  structure?: Structure
  meta: SupercellBuildMeta
}
