import { useState } from 'react'
import type { ReactNode } from 'react'
import {
  Activity,
  ArrowLeftRight,
  Atom,
  FileText,
  Grid3x3,
  Menu,
  Plus,
  Search,
  Settings,
  UserCircle,
} from 'lucide-react'

import type { ToolMode, WorkspaceFile } from './types'
import { FilePanel } from './components/FilePanel'
import { ToolPanel } from './components/ToolPanel'
import { DEMO_STRUCTURES } from './demoStructures'

const MOCK_FILES: WorkspaceFile[] = [
  {
    id: '1',
    name: 'benzen.in',
    kind: 'in',
    label: 'Benzen',
    pdbText: DEMO_STRUCTURES.benzene,
    initialOpenSections: { table: false, parameter: false },
  },
  {
    id: '2',
    name: 'h2o.in',
    kind: 'in',
    label: 'H2O',
    pdbText: DEMO_STRUCTURES.h2o,
    initialOpenSections: { table: false, parameter: true },
  },
  {
    id: '3',
    name: 'phenol.in',
    kind: 'in',
    label: 'Phenol',
    pdbText: DEMO_STRUCTURES.phenol,
    initialOpenSections: { table: false, parameter: false },
  },
]

const TOOL_NAV: Array<{ id: ToolMode; label: string; icon: ReactNode }> =
  [
    {
      id: 'transfer',
      label: 'Transfer',
      icon: <ArrowLeftRight />,
    },
    {
      id: 'supercell',
      label: 'Supercell',
      icon: <Grid3x3 />,
    },
    {
      id: 'vibration',
      label: 'Vibrations',
      icon: <Activity />,
    },
  ]

const FILES_BY_ID = new Map(MOCK_FILES.map((file) => [file.id, file]))

export default function EditorV2Page() {
  const [openFileIds, setOpenFileIds] = useState<string[]>(['1', '2', '3'])
  const [activeTool, setActiveTool] = useState<ToolMode | null>(null)

  const openFile = (id: string) => {
    if (openFileIds.includes(id)) {
      return
    }
    setOpenFileIds((prev) => [...prev, id])
  }

  const closeFile = (id: string) => {
    setOpenFileIds((prev) => prev.filter((fileId) => fileId !== id))
  }

  const toggleTool = (mode: ToolMode) => {
    setActiveTool((prev) => (prev === mode ? null : mode))
  }

  return (
    <div className="flex h-screen w-full overflow-hidden bg-background font-sans text-foreground">
      <aside className="z-20 flex w-16 flex-shrink-0 flex-col items-center gap-6 border-r border-border bg-slate-50 py-4">
        <button
          type="button"
          className="rounded-md p-2 transition-colors hover:bg-slate-200"
        >
          <Menu className="h-6 w-6 text-slate-700" />
        </button>

        <div className="mt-4 flex w-full flex-col items-center gap-6">
          {TOOL_NAV.map((item) => (
            <NavItem
              key={item.id}
              icon={item.icon}
              label={item.label}
              onClick={() => toggleTool(item.id)}
              isActive={activeTool === item.id}
            />
          ))}
        </div>

        <div className="mt-auto">
          <NavItem icon={<Settings />} label="Settings" />
        </div>
      </aside>

      <div className="flex min-w-0 flex-1 flex-col">
        <header className="z-10 flex h-16 items-center justify-between border-b border-border bg-white px-6 shadow-sm">
          <div className="flex items-center gap-3">
            <div className="flex h-8 w-8 items-center justify-center rounded-lg bg-blue-600 text-white">
              <Atom className="h-5 w-5" />
            </div>
            <span className="text-lg font-bold tracking-tight">QESpresso</span>
          </div>

          <div className="mx-8 flex-1">
            <div className="group relative max-w-xl">
              <Search className="absolute left-3 top-1/2 h-4 w-4 -translate-y-1/2 text-muted-foreground group-focus-within:text-blue-500" />
              <input
                type="text"
                placeholder="Search tools..."
                className="w-full rounded-md border-none bg-slate-100 py-2 pl-10 pr-4 text-sm text-slate-800 outline-none transition-all placeholder:text-muted-foreground focus:bg-white focus:ring-2 focus:ring-blue-500/20"
              />
            </div>
          </div>

          <button
            type="button"
            className="rounded-full p-1.5 transition-colors hover:bg-slate-100"
          >
            <UserCircle className="h-8 w-8 text-slate-600" />
          </button>
        </header>

        <div className="flex flex-1 overflow-hidden">
          <div className="flex w-64 flex-col border-r border-border bg-slate-50/50">
            <div className="border-b border-border bg-white p-4">
              <h2 className="text-sm font-semibold text-slate-900">Structures</h2>
            </div>

            <div className="flex-1 space-y-1 overflow-y-auto p-3">
              {MOCK_FILES.map((file) => {
                const isOpen = openFileIds.includes(file.id)
                return (
                  <button
                    type="button"
                    key={file.id}
                    onClick={() => (isOpen ? null : openFile(file.id))}
                    className={`flex w-full items-center gap-2 rounded-md border px-3 py-2 text-left text-sm transition-all ${
                      isOpen
                        ? 'border-blue-200 bg-blue-50 shadow-sm'
                        : 'border-slate-200 bg-white opacity-70 hover:border-blue-300 hover:opacity-100 hover:shadow-sm'
                    }`}
                  >
                    <FileText
                      className={`h-4 w-4 ${
                        isOpen ? 'text-blue-600' : 'text-slate-400'
                      }`}
                    />
                    <span
                      className={`truncate font-medium ${
                        isOpen ? 'text-blue-900' : 'text-slate-600'
                      }`}
                    >
                      {file.name}
                    </span>
                    {!isOpen ? (
                      <span className="ml-auto text-[10px] text-slate-400">
                        Hidden
                      </span>
                    ) : null}
                  </button>
                )
              })}

              <div className="m-2 rounded-lg border-2 border-dashed border-slate-200 bg-slate-50/50 px-4 py-8 text-center">
                <p className="text-xs text-muted-foreground">
                  Drag files here to import
                </p>
              </div>
            </div>

            <div className="border-t border-border bg-white p-4">
              <button className="flex w-full items-center justify-center gap-2 rounded-md bg-slate-900 py-2 text-white shadow-sm transition-colors hover:bg-slate-800">
                <Plus className="h-4 w-4" />
                <span className="text-sm font-medium">Import File</span>
              </button>
            </div>
          </div>

          <main className="flex-1 overflow-x-auto overflow-y-hidden bg-slate-100/50">
            <div className="flex h-full min-w-max items-start gap-6 p-6">
              {openFileIds.map((id) => {
                const fileData = FILES_BY_ID.get(id)
                if (!fileData) {
                  return null
                }
                return (
                  <div
                    key={id}
                    className="h-full overflow-hidden rounded-xl border border-slate-200 bg-white shadow-sm transition-shadow duration-300 hover:shadow-md"
                  >
                    <FilePanel data={fileData} onClose={() => closeFile(id)} />
                  </div>
                )
              })}

              {activeTool ? (
                <div className="h-full overflow-hidden rounded-xl border border-blue-100 bg-white shadow-lg ring-1 ring-blue-100">
                  <ToolPanel mode={activeTool} onClose={() => setActiveTool(null)} />
                </div>
              ) : null}

              <div className="h-full w-12 opacity-0" />
            </div>
          </main>
        </div>
      </div>
    </div>
  )
}

interface NavItemProps {
  icon: ReactNode
  label: string
  onClick?: () => void
  isActive?: boolean
}

function NavItem({ icon, label, onClick, isActive }: NavItemProps) {
  return (
    <button
      type="button"
      onClick={onClick}
      className="group flex w-full flex-col items-center gap-1 px-2"
    >
      <div
        className={`rounded-xl p-2 transition-all duration-200 ${
          isActive
            ? 'bg-blue-600 text-white shadow-md'
            : 'text-slate-500 hover:bg-blue-50 hover:text-blue-600'
        }`}
      >
        {icon}
      </div>
      <span
        className={`text-center text-[10px] font-medium leading-none ${
          isActive
            ? 'font-bold text-blue-700'
            : 'text-slate-500 group-hover:text-blue-700'
        }`}
      >
        {label}
      </span>
    </button>
  )
}
