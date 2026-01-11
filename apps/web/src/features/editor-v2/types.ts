import type { Structure } from '@/lib/types'

export type WorkspaceFile = {
  id: string
  name: string
  kind: 'in' | 'out'
  label: string
  pdbText?: string
  bcifUrl?: string
  structureId?: string
  structure?: Structure
  parseSource?: string
  initialOpenSections: {
    table: boolean
    parameter: boolean
  }
}

export type ToolMode = 'transfer' | 'supercell' | 'vibration'
