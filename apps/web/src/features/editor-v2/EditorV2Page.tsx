import { useCallback, useEffect, useMemo, useRef, useState } from 'react'
import type { ReactNode } from 'react'
import type { DockviewApi, DockviewReadyEvent, IDockviewPanelProps } from 'dockview-react'
import { DockviewReact } from 'dockview-react'
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
import 'dockview/dist/styles/dockview.css'
import { parseQeInput } from '@/lib/api'
import { atomsToPdb } from '@/lib/pdb'

const INITIAL_FILES: WorkspaceFile[] = []

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

const FILE_PANEL_PREFIX = 'file'
const TOOL_PANEL_PREFIX = 'tool'
const filePanelId = (id: string) => `${FILE_PANEL_PREFIX}-${id}`
const toolPanelId = (id: ToolMode) => `${TOOL_PANEL_PREFIX}-${id}`

export default function EditorV2Page() {
  const [files, setFiles] = useState<WorkspaceFile[]>(() => [...INITIAL_FILES])
  const [importError, setImportError] = useState<string | null>(null)
  const [isImporting, setIsImporting] = useState(false)
  const [activeTool, setActiveTool] = useState<ToolMode | null>(null)
  const [activeFileId, setActiveFileId] = useState<string | null>(
    INITIAL_FILES[0]?.id ?? null,
  )
  const [pendingOpenFileId, setPendingOpenFileId] = useState<string | null>(
    null,
  )
  const dockviewApiRef = useRef<DockviewApi | null>(null)
  const disposablesRef = useRef<Array<{ dispose: () => void }>>([])
  const dockviewContainerRef = useRef<HTMLDivElement | null>(null)
  const fileInputRef = useRef<HTMLInputElement | null>(null)

  const filesById = useMemo(
    () => new Map(files.map((file) => [file.id, file])),
    [files],
  )

  const openFile = useCallback((id: string) => {
    const api = dockviewApiRef.current
    const file = filesById.get(id)
    if (!api || !file) {
      return
    }
    const panelId = filePanelId(id)
    const existing = api.getPanel(panelId)
    if (existing) {
      existing.api.setActive()
      return
    }
    api.addPanel({
      id: panelId,
      title: file.name,
      component: 'structure',
      params: { fileId: id },
    })
  }, [filesById])

  const openTool = useCallback((mode: ToolMode) => {
    const api = dockviewApiRef.current
    if (!api) {
      return
    }
    const panelId = toolPanelId(mode)
    const existing = api.getPanel(panelId)
    if (existing) {
      existing.api.setActive()
      return
    }
    api.addPanel({
      id: panelId,
      title: TOOL_NAV.find((tool) => tool.id === mode)?.label ?? mode,
      component: 'tool',
      params: { mode },
    })
  }, [])

  const handleReady = useCallback((event: DockviewReadyEvent) => {
    dockviewApiRef.current = event.api

    disposablesRef.current = [
      event.api.onDidActivePanelChange((panel) => {
        if (!panel) {
          setActiveFileId(null)
          setActiveTool(null)
          return
        }
        if (panel.id.startsWith(`${FILE_PANEL_PREFIX}-`)) {
          setActiveFileId(panel.id.replace(`${FILE_PANEL_PREFIX}-`, ''))
          setActiveTool(null)
          return
        }
        if (panel.id.startsWith(`${TOOL_PANEL_PREFIX}-`)) {
          const mode = panel.id.replace(`${TOOL_PANEL_PREFIX}-`, '') as ToolMode
          setActiveTool(mode)
          return
        }
        setActiveTool(null)
      }),
    ]

    requestAnimationFrame(() => {
      const container = dockviewContainerRef.current
      if (container) {
        event.api.layout(container.clientWidth, container.clientHeight, true)
      }
    })
  }, [])

  const dockviewComponents = useMemo(
    () => ({
      structure: ({
        params,
      }: IDockviewPanelProps<{ fileId: string }>) => {
        const file = params?.fileId ? filesById.get(params.fileId) : null
        if (!file) {
          return (
            <div className="flex h-full items-center justify-center text-sm text-muted-foreground">
              Missing structure data
            </div>
          )
        }
        return (
          <FilePanel
            data={file}
            showHeader={false}
            className="h-full w-full border-none p-3"
          />
        )
      },
      tool: ({
        params,
        api,
      }: IDockviewPanelProps<{ mode: ToolMode }>) => {
        if (!params?.mode) {
          return null
        }
        return (
          <ToolPanel
            mode={params.mode}
            onClose={() => api.close()}
            variant="dock"
            showClose={false}
            className="h-full w-full border-none"
          />
        )
      },
      history: () => <HistoryPanel />,
    }),
    [filesById],
  )

  const handleImportFile = useCallback(
    async (event: React.ChangeEvent<HTMLInputElement>) => {
      const file = event.target.files?.[0]
      if (!file) {
        return
      }
      setImportError(null)
      setIsImporting(true)
      try {
        const content = await file.text()
        const structure = await parseQeInput(content)
        const baseName = file.name.replace(/\.[^/.]+$/, '') || file.name
        const id = `import-${Date.now()}-${Math.random()
          .toString(36)
          .slice(2, 8)}`
        const pdbText =
          structure.atoms.length > 0 ? atomsToPdb(structure.atoms) : undefined
        const nextFile: WorkspaceFile = {
          id,
          name: file.name,
          kind: 'in',
          label: baseName,
          pdbText,
          initialOpenSections: { table: false, parameter: true },
        }
        setFiles((prev) => [...prev, nextFile])
        setActiveFileId(id)
        setPendingOpenFileId(id)
      } catch (err) {
        setImportError(
          err instanceof Error ? err.message : 'Failed to import file.',
        )
      } finally {
        setIsImporting(false)
        event.target.value = ''
      }
    },
    [],
  )

  useEffect(() => {
    if (!pendingOpenFileId) {
      return
    }
    if (!filesById.has(pendingOpenFileId)) {
      return
    }
    openFile(pendingOpenFileId)
    setPendingOpenFileId(null)
  }, [filesById, openFile, pendingOpenFileId])

  useEffect(() => {
    const container = dockviewContainerRef.current
    if (!container) {
      return () => {
        disposablesRef.current.forEach((disposable) => disposable.dispose())
        disposablesRef.current = []
      }
    }

    const observer = new ResizeObserver((entries) => {
      const entry = entries[0]
      if (!entry) {
        return
      }
      const api = dockviewApiRef.current
      if (!api) {
        return
      }
      api.layout(entry.contentRect.width, entry.contentRect.height, true)
    })
    observer.observe(container)

    return () => {
      observer.disconnect()
      disposablesRef.current.forEach((disposable) => disposable.dispose())
      disposablesRef.current = []
    }
  }, [])

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
              onClick={() => {
                setActiveTool(item.id)
                openTool(item.id)
              }}
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
              {files.map((file) => {
                const isActive = activeFileId === file.id
                return (
                  <button
                    type="button"
                    key={file.id}
                    onClick={() => openFile(file.id)}
                    className={`flex w-full items-center gap-2 rounded-md border px-3 py-2 text-left text-sm transition-all ${
                      isActive
                        ? 'border-blue-200 bg-blue-50 shadow-sm'
                        : 'border-slate-200 bg-white opacity-70 hover:border-blue-300 hover:opacity-100 hover:shadow-sm'
                    }`}
                  >
                    <FileText
                      className={`h-4 w-4 ${
                        isActive ? 'text-blue-600' : 'text-slate-400'
                      }`}
                    />
                    <span
                      className={`truncate font-medium ${
                        isActive ? 'text-blue-900' : 'text-slate-600'
                      }`}
                    >
                      {file.name}
                    </span>
                    {isActive ? (
                      <span className="ml-auto text-[10px] text-blue-500">
                        Active
                      </span>
                    ) : null}
                  </button>
                )
              })}

              {files.length === 0 ? (
                <div className="rounded-md border border-dashed border-slate-200 bg-white px-3 py-4 text-xs text-slate-500">
                  No files imported yet.
                </div>
              ) : null}

              <div className="m-2 rounded-lg border-2 border-dashed border-slate-200 bg-slate-50/50 px-4 py-8 text-center">
                <p className="text-xs text-muted-foreground">
                  Drag files here to import
                </p>
              </div>
            </div>

            <div className="border-t border-border bg-white p-4">
              {importError ? (
                <p className="mb-2 rounded-md bg-red-50 px-3 py-2 text-xs text-red-600">
                  {importError}
                </p>
              ) : null}
              <button
                type="button"
                onClick={() => fileInputRef.current?.click()}
                disabled={isImporting}
                className="flex w-full items-center justify-center gap-2 rounded-md bg-slate-900 py-2 text-white shadow-sm transition-colors hover:bg-slate-800 disabled:cursor-not-allowed disabled:opacity-60"
              >
                <Plus className="h-4 w-4" />
                <span className="text-sm font-medium">
                  {isImporting ? 'Importing…' : 'Import File'}
                </span>
              </button>
              <input
                ref={fileInputRef}
                type="file"
                accept=".in"
                onChange={handleImportFile}
                className="hidden"
              />
            </div>
          </div>

          <main className="flex-1 overflow-hidden bg-slate-100/50 p-4 min-h-0">
            <div
              ref={dockviewContainerRef}
              className="relative h-full w-full min-h-0"
            >
              <DockviewReact
                components={dockviewComponents}
                onReady={handleReady}
                className="dockview-theme-light h-full w-full"
              />
              {files.length === 0 ? (
                <div className="pointer-events-none absolute inset-0 flex flex-col items-center justify-center gap-3 text-sm text-slate-500">
                  <div className="rounded-full bg-white px-4 py-2 shadow-sm">
                    Import a .in file to start.
                  </div>
                  <span className="text-xs text-slate-400">
                    Your workspace begins empty.
                  </span>
                </div>
              ) : null}
            </div>
          </main>
        </div>
      </div>
    </div>
  )
}

function HistoryPanel() {
  return (
    <div className="flex h-full flex-col gap-4 bg-white p-4 text-sm text-slate-700">
      <div className="flex items-center justify-between">
        <div>
          <p className="text-xs uppercase tracking-[0.2em] text-slate-400">
            Lineage
          </p>
          <p className="text-base font-semibold text-slate-800">
            Derived Structures
          </p>
        </div>
        <span className="rounded-full bg-slate-100 px-2 py-1 text-xs text-slate-500">
          3 items
        </span>
      </div>
      <div className="flex-1 space-y-3 overflow-y-auto">
        {['benzen.in → transfer', 'h2o.in → supercell', 'phenol.in → draft'].map(
          (item, index) => (
            <div
              key={`${item}-${index}`}
              className="rounded-lg border border-slate-200 bg-slate-50 px-3 py-2"
            >
              <p className="text-xs font-medium text-slate-600">Step {index + 1}</p>
              <p className="text-sm text-slate-800">{item}</p>
              <p className="text-[10px] text-slate-400">Pending review</p>
            </div>
          ),
        )}
      </div>
      <div className="rounded-lg border border-dashed border-slate-200 bg-slate-50 px-3 py-4 text-center text-xs text-slate-400">
        Drag a panel to link lineage
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
