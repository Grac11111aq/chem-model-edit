export type Atom = {
  symbol: string
  x: number
  y: number
  z: number
}

export type Structure = {
  atoms: Atom[]
}

export type Vector3 = {
  x: number
  y: number
  z: number
}

export type SupercellMeta = {
  na: number
  nb: number
  layers: number
}
