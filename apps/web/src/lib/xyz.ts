import type { Atom } from './types'

const LINE_RE =
  /^\s*(?<sym>[A-Za-z]{1,3})\s+(?<x>-?\d+\.?\d*(?:[eE][+-]?\d+)?)\s+(?<y>-?\d+\.?\d*(?:[eE][+-]?\d+)?)\s+(?<z>-?\d+\.?\d*(?:[eE][+-]?\d+)?)/

export function parseXyzBlock(text: string): Atom[] {
  const atoms: Atom[] = []
  const lines = text.split(/\r?\n/).map((line) => line.trim())
  for (const line of lines) {
    if (!line) continue
    const match = line.match(LINE_RE)
    if (!match?.groups) continue
    atoms.push({
      symbol: match.groups.sym,
      x: Number(match.groups.x),
      y: Number(match.groups.y),
      z: Number(match.groups.z),
    })
  }
  return atoms
}

export function atomsToXyz(atoms: Atom[]): string {
  return atoms
    .map(
      (atom) =>
        `${atom.symbol} ${atom.x.toFixed(6)} ${atom.y.toFixed(6)} ${atom.z.toFixed(6)}`,
    )
    .join('\n')
}
