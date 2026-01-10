import { Cuboid, Minus, X } from 'lucide-react'

import type { WorkspaceFile } from '../types'
import { CollapsibleSection } from './CollapsibleSection'
import MolstarViewer from '@/components/molstar/MolstarViewer'

interface FilePanelProps {
  data: WorkspaceFile
  onClose: () => void
}

export function FilePanel({ data, onClose }: FilePanelProps) {
  return (
    <div className="flex h-full w-80 flex-shrink-0 flex-col gap-4 border-r border-border bg-background p-4 animate-in fade-in slide-in-from-left-4 duration-300">
      <div className="flex items-center justify-between">
        <h2
          className="flex-1 truncate text-lg font-semibold tracking-tight"
          title={data.name}
        >
          {data.name}
        </h2>
        <div className="flex items-center gap-1">
          <button
            type="button"
            onClick={onClose}
            className="rounded p-1 text-slate-400 transition-colors hover:bg-slate-100 hover:text-slate-600"
            title="Minimize"
          >
            <Minus className="h-4 w-4" />
          </button>
          <button
            type="button"
            onClick={onClose}
            className="rounded p-1 text-slate-400 transition-colors hover:bg-red-50 hover:text-red-500"
            title="Close"
          >
            <X className="h-4 w-4" />
          </button>
        </div>
      </div>

      <div className="group relative flex w-full flex-col overflow-hidden rounded-lg border border-slate-200 bg-white">
        <div className="absolute inset-0 bg-grid-slate-200/50 [mask-image:linear-gradient(0deg,white,rgba(255,255,255,0.7))]" />
        <div className="relative z-10 flex h-48 w-full flex-col">
          {data.pdbText ? (
            <MolstarViewer pdbText={data.pdbText} />
          ) : (
            <div className="flex h-full w-full flex-col items-center justify-center px-4 text-center text-muted-foreground">
              <Cuboid className="mb-2 h-12 w-12 text-slate-300 transition-colors duration-300 group-hover:text-blue-400" />
              <span className="font-medium text-slate-600">{data.label}</span>
              <span className="text-xs text-slate-400">3D Structure View</span>
            </div>
          )}
        </div>
        <div className="flex items-center justify-between border-t border-slate-200 bg-white px-3 py-2 text-xs text-slate-500">
          <span className="font-medium">{data.label}</span>
          <span className="rounded-full bg-slate-100 px-2 py-0.5 text-[10px] uppercase tracking-wide">
            Mol*
          </span>
        </div>
      </div>

      <div className="flex flex-1 flex-col gap-2 overflow-y-auto pr-1">
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
              <div className="grid grid-cols-4 gap-2">
                <span>C</span> <span>1.45</span> <span>3.62</span>{' '}
                <span>4.25</span>
              </div>
              <div className="grid grid-cols-4 gap-2">
                <span>H</span> <span>0.82</span> <span>2.11</span>{' '}
                <span>3.10</span>
              </div>
              <div className="grid grid-cols-4 gap-2">
                <span>O</span> <span>2.15</span> <span>4.01</span>{' '}
                <span>1.22</span>
              </div>
            </div>
          </div>
        </CollapsibleSection>

        <CollapsibleSection
          title="Parameters"
          defaultOpen={data.initialOpenSections.parameter}
        >
          <div className="space-y-2">
            <p className="mb-1 text-xs text-muted-foreground">QE Option Params</p>
            <div className="grid grid-cols-2 gap-2 text-xs">
              <div className="flex flex-col gap-1">
                <label className="text-slate-500">ecutwfc</label>
                <input
                  type="text"
                  value="30.0"
                  className="rounded border border-slate-200 bg-slate-50 px-2 py-1"
                  readOnly
                />
              </div>
              <div className="flex flex-col gap-1">
                <label className="text-slate-500">mixing_beta</label>
                <input
                  type="text"
                  value="0.7"
                  className="rounded border border-slate-200 bg-slate-50 px-2 py-1"
                  readOnly
                />
              </div>
              <div className="flex flex-col gap-1">
                <label className="text-slate-500">conv_thr</label>
                <input
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
  )
}
