export type WorkspaceFile = {
  id: string
  name: string
  kind: 'in' | 'out'
  label: string
  pdbText?: string
  initialOpenSections: {
    table: boolean
    parameter: boolean
  }
}

export type ToolMode = 'transfer' | 'supercell' | 'vibration'
