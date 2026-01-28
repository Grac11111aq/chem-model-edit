import { getAuthToken } from './auth'
import type {
  AuthMe,
  AuthSession,
  Lattice,
  LatticeParams,
  Structure,
  SupercellBuildRequest,
  SupercellBuildResponse,
  SupercellMeta,
  ZPEJobRequest,
  ZPEJobResponse,
  ZPEJobStatus,
  ZPEParseResponse,
  ZPEQueueTargetList,
  ZPEResult,
} from './types'

const API_BASE = import.meta.env.VITE_API_BASE ?? 'http://localhost:8000'

type ApiError = {
  detail?: string
}

const withAuthHeaders = (headers: Record<string, string> = {}) => {
  const token = getAuthToken()
  if (!token) {
    return headers
  }
  return { ...headers, Authorization: `Bearer ${token}` }
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

async function handleTextResponse(response: Response): Promise<string> {
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
  return response.text()
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

export async function createStructureFromQe(content: string): Promise<{
  structure_id: string
  structure: Structure
  source: string
}> {
  const response = await fetch(`${API_BASE}/structures`, {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify({ content }),
  })
  return handleResponse<{
    structure_id: string
    structure: Structure
    source: string
  }>(response)
}

export async function getStructure(structureId: string): Promise<Structure> {
  const safeId = encodeURIComponent(structureId)
  const response = await fetch(`${API_BASE}/structures/${safeId}`)
  const data = await handleResponse<{ structure: Structure }>(response)
  return data.structure
}

export function structureViewUrl(
  structureId: string,
  params?: {
    format?: 'cif'
  },
): string {
  const format = params?.format ?? 'cif'
  const query = new URLSearchParams({
    format,
  })
  const safeId = encodeURIComponent(structureId)
  return `${API_BASE}/structures/${safeId}/view?${query.toString()}`
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

export async function buildSupercell(
  params: SupercellBuildRequest,
): Promise<SupercellBuildResponse> {
  const response = await fetch(`${API_BASE}/supercell/build`, {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(params),
  })
  return handleResponse<SupercellBuildResponse>(response)
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

export async function parseZpeInput(
  content: string,
  structureId?: string | null,
): Promise<ZPEParseResponse> {
  const payload: { content: string; structure_id?: string } = { content }
  if (structureId) {
    payload.structure_id = structureId
  }
  const response = await fetch(`${API_BASE}/calc/zpe/parse`, {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(payload),
  })
  return handleResponse<ZPEParseResponse>(response)
}

export async function registerAccount(params: {
  email: string
  password: string
}): Promise<AuthSession> {
  const response = await fetch(`${API_BASE}/auth/register`, {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(params),
  })
  return handleResponse<AuthSession>(response)
}

export async function loginAccount(params: {
  email: string
  password: string
}): Promise<AuthSession> {
  const response = await fetch(`${API_BASE}/auth/login`, {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(params),
  })
  return handleResponse<AuthSession>(response)
}

export async function fetchAuthMe(): Promise<AuthMe> {
  const response = await fetch(`${API_BASE}/auth/me`, {
    headers: withAuthHeaders(),
  })
  return handleResponse<AuthMe>(response)
}

export async function logoutAccount(): Promise<void> {
  const response = await fetch(`${API_BASE}/auth/logout`, {
    method: 'POST',
    headers: withAuthHeaders(),
  })
  await handleResponse<{ ok: boolean }>(response)
}

export async function createEnrollToken(params?: {
  ttlSeconds?: number
  label?: string
}): Promise<{
  token: string
  expires_at: string
  ttl_seconds: number
  label?: string | null
}> {
  const response = await fetch(`${API_BASE}/calc/zpe/compute/enroll-tokens`, {
    method: 'POST',
    headers: withAuthHeaders({ 'Content-Type': 'application/json' }),
    body: JSON.stringify({
      ttl_seconds: params?.ttlSeconds,
      label: params?.label,
    }),
  })
  return handleResponse(response)
}

export async function fetchQueueTargets(): Promise<ZPEQueueTargetList> {
  const response = await fetch(`${API_BASE}/calc/zpe/compute/targets`, {
    headers: withAuthHeaders(),
  })
  return handleResponse<ZPEQueueTargetList>(response)
}

export async function selectQueueTarget(
  targetId: string,
): Promise<{ active_target_id: string }> {
  const response = await fetch(`${API_BASE}/calc/zpe/compute/targets/select`, {
    method: 'POST',
    headers: withAuthHeaders({ 'Content-Type': 'application/json' }),
    body: JSON.stringify({ target_id: targetId }),
  })
  return handleResponse(response)
}

export async function createZpeJob(
  request: ZPEJobRequest,
): Promise<ZPEJobResponse> {
  const response = await fetch(`${API_BASE}/calc/zpe/jobs`, {
    method: 'POST',
    headers: withAuthHeaders({ 'Content-Type': 'application/json' }),
    body: JSON.stringify(request),
  })
  return handleResponse<ZPEJobResponse>(response)
}

export async function fetchZpeStatus(jobId: string): Promise<ZPEJobStatus> {
  const safeId = encodeURIComponent(jobId)
  const response = await fetch(`${API_BASE}/calc/zpe/jobs/${safeId}`, {
    headers: withAuthHeaders(),
  })
  return handleResponse<ZPEJobStatus>(response)
}

export async function fetchZpeResult(jobId: string): Promise<ZPEResult> {
  const safeId = encodeURIComponent(jobId)
  const response = await fetch(`${API_BASE}/calc/zpe/jobs/${safeId}/result`, {
    headers: withAuthHeaders(),
  })
  const data = await handleResponse<{ result: ZPEResult }>(response)
  return data.result
}

export async function downloadZpeFile(
  jobId: string,
  kind: 'summary' | 'freqs',
): Promise<string> {
  const safeId = encodeURIComponent(jobId)
  const response = await fetch(
    `${API_BASE}/calc/zpe/jobs/${safeId}/files?kind=${kind}`,
    { headers: withAuthHeaders() },
  )
  return handleTextResponse(response)
}
