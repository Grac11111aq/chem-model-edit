import { useEffect, useRef, useState } from 'react'
import { Cuboid, Minus, X } from 'lucide-react'

import { CollapsibleSection } from './CollapsibleSection'

import type { WorkspaceFile } from '../types'
import type { Structure } from '@/lib/types'

import MolstarViewer from '@/components/molstar/MolstarViewer'
import { getStructure } from '@/lib/api'
import { cn } from '@/lib/utils'

interface FilePanelProps {
  data: WorkspaceFile
  fileId: string
  onStructureLoaded?: (fileId: string, structure: Structure) => void
  onClose?: () => void
  onMinimize?: () => void
  showHeader?: boolean
  className?: string
}

export function FilePanel({
  data,
  fileId,
  onStructureLoaded,
  onClose,
  onMinimize,
  showHeader = true,
  className,
}: FilePanelProps) {
  const [viewerError, setViewerError] = useState<string | null>(null)
  const [structure, setStructure] = useState<Structure | null>(
    data.structure ?? null,
  )
  const [tableError, setTableError] = useState<string | null>(null)
  const [isTableLoading, setIsTableLoading] = useState(false)
  const onStructureLoadedRef = useRef<FilePanelProps['onStructureLoaded']>(
    onStructureLoaded,
  )

  useEffect(() => {
    setViewerError(null)
  }, [data.cifUrl, data.pdbText])

  useEffect(() => {
    onStructureLoadedRef.current = onStructureLoaded
  }, [onStructureLoaded])

  useEffect(() => {
    setStructure(data.structure ?? null)
  }, [data.structure, data.structureId])

  useEffect(() => {
    if (data.structure) {
      setIsTableLoading(false)
      setTableError(null)
      return
    }
    if (!data.structureId) {
      setIsTableLoading(false)
      setTableError(null)
      return
    }
    let cancelled = false
    setIsTableLoading(true)
    setTableError(null)
    getStructure(data.structureId)
      .then((nextStructure) => {
        if (cancelled) return
        setStructure(nextStructure)
        setIsTableLoading(false)
        onStructureLoadedRef.current?.(fileId, nextStructure)
      })
      .catch((err) => {
        if (cancelled) return
        setTableError(err instanceof Error ? err.message : String(err))
        setIsTableLoading(false)
      })
    return () => {
      cancelled = true
    }
  }, [data.structure, data.structureId])

  const atoms = structure?.atoms ?? []

  return (
    <div
      className={cn(
        'flex h-full w-80 flex-shrink-0 flex-col gap-4 border-r border-border bg-background p-4 animate-in fade-in slide-in-from-left-4 duration-300',
        className,
      )}
    >
      {showHeader ? (
        <div className="flex items-center justify-between">
          <h2
            className="flex-1 truncate text-lg font-semibold tracking-tight"
            title={data.name}
          >
            {data.name}
          </h2>
          {onClose || onMinimize ? (
            <div className="flex items-center gap-1">
              {onMinimize ? (
                <button
                  type="button"
                  onClick={onMinimize}
                  className="rounded p-1 text-slate-400 transition-colors hover:bg-slate-100 hover:text-slate-600"
                  title="Minimize"
                  aria-label="Minimize file panel"
                >
                  <Minus className="h-4 w-4" aria-hidden="true" />
                </button>
              ) : null}
              {onClose ? (
                <button
                  type="button"
                  onClick={onClose}
                  className="rounded p-1 text-slate-400 transition-colors hover:bg-red-50 hover:text-red-500"
                  title="Close"
                  aria-label="Close file panel"
                >
                  <X className="h-4 w-4" aria-hidden="true" />
                </button>
              ) : null}
            </div>
          ) : null}
        </div>
      ) : null}

      <div className="flex min-h-0 flex-1 flex-col gap-4">
        <div className="group relative flex h-1/2 min-h-[220px] w-full flex-col overflow-hidden rounded-lg border border-border bg-card">
          <div className="absolute inset-0 bg-grid-slate-200/50 [mask-image:linear-gradient(0deg,white,rgba(255,255,255,0.7))]" />
          <div className="relative z-10 flex w-full flex-1 flex-col">
            {data.cifUrl ? (
              <MolstarViewer
                cifUrl={data.cifUrl}
                onError={setViewerError}
                onLoad={() => setViewerError(null)}
              />
            ) : data.pdbText ? (
              <MolstarViewer
                pdbText={data.pdbText}
                onError={setViewerError}
                onLoad={() => setViewerError(null)}
              />
            ) : (
              <div className="flex h-full w-full flex-col items-center justify-center px-4 text-center text-muted-foreground">
                <Cuboid className="mb-2 h-12 w-12 text-slate-300 transition-colors duration-300 group-hover:text-blue-400" />
                <span className="font-medium text-slate-600">{data.label}</span>
                <span className="text-xs text-slate-400">
                  3D Structure View
                </span>
              </div>
            )}
          </div>
          {viewerError ? (
            <div className="absolute inset-0 z-20 flex items-center justify-center bg-white/80 px-4 text-center">
              <div className="rounded-md border border-red-200 bg-red-50 px-3 py-2 text-sm text-red-700 shadow-sm">
                <p className="font-semibold">Viewer failed to load</p>
                <p className="mt-1 text-xs text-red-600">{viewerError}</p>
              </div>
            </div>
          ) : null}
        </div>

        <div className="flex min-h-0 flex-1 flex-col gap-2 overflow-y-auto pr-1">
          <CollapsibleSection
            title="Table"
            defaultOpen={data.initialOpenSections.table}
          >
            <div className="space-y-2">
              <p className="mb-1 text-xs text-muted-foreground">
                Atomic Positions (Angstrom)
              </p>
              <div className="rounded border border-slate-100 bg-slate-50 p-2 text-xs font-mono">
                <div className="mb-1 grid grid-cols-4 gap-2 border-b border-slate-200 pb-1 font-bold text-slate-600">
                  <span>El</span>
                  <span>X</span>
                  <span>Y</span>
                  <span>Z</span>
                </div>
                {isTableLoading ? (
                  <div className="py-2 text-slate-400">Loading...</div>
                ) : tableError ? (
                  <div className="py-2 text-red-500">
                    Failed to load structure.
                  </div>
                ) : atoms.length === 0 ? (
                  <div className="py-2 text-slate-400">No atoms.</div>
                ) : (
                  atoms.map((atom, index) => (
                    <div key={`${atom.symbol}-${index}`} className="grid grid-cols-4 gap-2">
                      <span>{atom.symbol}</span>
                      <span>{atom.x.toFixed(4)}</span>
                      <span>{atom.y.toFixed(4)}</span>
                      <span>{atom.z.toFixed(4)}</span>
                    </div>
                  ))
                )}
              </div>
            </div>
          </CollapsibleSection>

          <CollapsibleSection
            title="Parameters"
            defaultOpen={data.initialOpenSections.parameter}
          >
            <div className="space-y-2">
              <p className="mb-1 text-xs text-muted-foreground">
                QE Option Params
              </p>
              <div className="grid grid-cols-2 gap-2 text-xs">
                <div className="flex flex-col gap-1">
                  <label htmlFor="param-ecutwfc" className="text-slate-500">
                    ecutwfc
                  </label>
                  <input
                    id="param-ecutwfc"
                    type="text"
                    value="30.0"
                    className="rounded border border-slate-200 bg-slate-50 px-2 py-1"
                    readOnly
                  />
                </div>
                <div className="flex flex-col gap-1">
                  <label htmlFor="param-mixing-beta" className="text-slate-500">
                    mixing_beta
                  </label>
                  <input
                    id="param-mixing-beta"
                    type="text"
                    value="0.7"
                    className="rounded border border-slate-200 bg-slate-50 px-2 py-1"
                    readOnly
                  />
                </div>
                <div className="flex flex-col gap-1">
                  <label htmlFor="param-conv-thr" className="text-slate-500">
                    conv_thr
                  </label>
                  <input
                    id="param-conv-thr"
                    type="text"
                    value="1.0d-8"
                    className="rounded border border-slate-200 bg-slate-50 px-2 py-1"
                    readOnly
                  />
                </div>
              </div>
            </div>
          </CollapsibleSection>
        </div>
      </div>
    </div>
  )
}
