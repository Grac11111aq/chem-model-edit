import type { Lattice, LatticeParams, Structure, SupercellMeta } from './types'

const API_BASE = import.meta.env.VITE_API_BASE ?? 'http://localhost:8000'

type ApiError = {
  detail?: string
}

async function handleResponse<T>(response: Response): Promise<T> {
  if (!response.ok) {
    let message = response.statusText
    try {
      const data = (await response.json()) as ApiError
      if (data.detail) {
        message = data.detail
      }
    } catch (_err) {
      // ignore JSON parsing errors
    }
    throw new Error(message)
  }
  return (await response.json()) as T
}

export async function parseQeInput(content: string): Promise<Structure> {
  const response = await fetch(`${API_BASE}/parse`, {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify({ content }),
  })
  const data = await handleResponse<{ structure: Structure }>(response)
  return data.structure
}

export async function exportQeInput(structure: Structure): Promise<string> {
  const response = await fetch(`${API_BASE}/export`, {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify({ structure }),
  })
  const data = await handleResponse<{ content: string }>(response)
  return data.content
}

export async function deltaTransplant(params: {
  smallIn: string
  smallOut: string
  largeIn: string
}): Promise<string> {
  const response = await fetch(`${API_BASE}/transplant/delta`, {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify({
      small_in: params.smallIn,
      small_out: params.smallOut,
      large_in: params.largeIn,
    }),
  })
  const data = await handleResponse<{ content: string }>(response)
  return data.content
}

export async function generateSupercell(params: {
  structureA: Structure
  structureB: Structure
  sequence: string
  lattice: Lattice
}): Promise<{ structure: Structure; meta: SupercellMeta }> {
  const response = await fetch(`${API_BASE}/supercell`, {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(params),
  })
  return handleResponse<{ structure: Structure; meta: SupercellMeta }>(response)
}

export async function generateTiledSupercell(params: {
  structureA: Structure
  structureB: Structure
  pattern: Array<Array<string>>
  lattice: Lattice
  checkOverlap: boolean
  overlapTolerance?: number
}): Promise<{ structure: Structure; meta: SupercellMeta }> {
  const response = await fetch(`${API_BASE}/supercell/tiled`, {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(params),
  })
  return handleResponse<{ structure: Structure; meta: SupercellMeta }>(response)
}

export async function latticeVectorsToParams(params: {
  lattice: Lattice
  unit: string
}): Promise<{ lattice: Lattice; params: LatticeParams; unit: string }> {
  const response = await fetch(`${API_BASE}/lattice/vectors-to-params`, {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(params),
  })
  return handleResponse<{
    lattice: Lattice
    params: LatticeParams
    unit: string
  }>(response)
}

export async function latticeParamsToVectors(params: {
  params: LatticeParams
  unit: string
}): Promise<{ lattice: Lattice; params: LatticeParams; unit: string }> {
  const response = await fetch(`${API_BASE}/lattice/params-to-vectors`, {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(params),
  })
  return handleResponse<{
    lattice: Lattice
    params: LatticeParams
    unit: string
  }>(response)
}
