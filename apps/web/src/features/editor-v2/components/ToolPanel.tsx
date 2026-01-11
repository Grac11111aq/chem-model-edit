import type { ReactNode } from 'react'
import { Activity, Grid3x3, Layers, MousePointerClick, Play, X } from 'lucide-react'

import type { ToolMode } from '../types'
import { CollapsibleSection } from './CollapsibleSection'
import { cn } from '@/lib/utils'

interface ToolPanelProps {
  mode: ToolMode
  onClose?: () => void
  variant?: 'stack' | 'dock'
  showHeader?: boolean
  showClose?: boolean
  className?: string
}

const toolTitles: Record<ToolMode, string> = {
  transfer: 'Structure Transfer',
  supercell: 'Supercell Builder',
  vibration: 'Vibrations (Preview)',
}

const toolIcons: Record<ToolMode, ReactNode> = {
  transfer: <Layers className="mb-2 h-16 w-16 text-emerald-300" />,
  supercell: <Grid3x3 className="mb-2 h-16 w-16 text-indigo-300" />,
  vibration: <Play className="mb-2 h-16 w-16 text-rose-300" />,
}

export function ToolPanel({
  mode,
  onClose,
  variant = 'stack',
  showHeader = true,
  showClose = true,
  className,
}: ToolPanelProps) {
  return (
    <div
      className={cn(
        'flex h-full flex-col gap-4 bg-slate-50/80 p-4',
        variant === 'stack'
          ? 'w-[60rem] flex-shrink-0 border-r border-border border-l-4 border-l-blue-500/20 animate-in fade-in slide-in-from-right-4 duration-300'
          : 'w-full',
        className,
      )}
    >
      {showHeader ? (
        <div className="flex shrink-0 items-center justify-between">
          <h2 className="text-lg font-semibold tracking-tight text-slate-800">
            {toolTitles[mode]}
          </h2>
          {showClose && onClose ? (
            <button
              type="button"
              onClick={onClose}
              className="rounded p-1 text-slate-400 transition-colors hover:bg-slate-200 hover:text-slate-600"
            >
              <X className="h-4 w-4" />
            </button>
          ) : null}
        </div>
      ) : null}

      <div className="flex flex-1 gap-4 overflow-hidden">
        <div className="flex w-[26rem] flex-shrink-0 flex-col gap-4 overflow-y-auto pr-1">
          <div className="relative flex aspect-square w-full flex-col items-center justify-center overflow-hidden rounded-lg border-2 border-dashed border-slate-200 bg-white text-muted-foreground shadow-sm">
            <div className="absolute inset-0 bg-grid-slate-100/50" />
            <div className="relative z-10 flex flex-col items-center">
              {toolIcons[mode]}
              <span className="font-medium text-slate-500">Action Preview</span>
              {mode === 'vibration' ? (
                <span className="mt-1 text-xs text-slate-400">Preview only</span>
              ) : null}
            </div>
          </div>

          <div className="flex flex-col gap-2 opacity-50 grayscale-[0.5]">
            <CollapsibleSection title="Table" defaultOpen={false}>
              <div className="h-20 rounded bg-slate-100" />
            </CollapsibleSection>
            <CollapsibleSection title="Parameters" defaultOpen={false}>
              <div className="h-20 rounded bg-slate-100" />
            </CollapsibleSection>
          </div>
        </div>

        <div className="flex min-w-0 flex-1 flex-col overflow-y-auto">
          <div className="flex h-full flex-col overflow-hidden rounded-md border border-blue-200 bg-white shadow-sm">
            <div className="flex shrink-0 items-center gap-2 border-b border-blue-100 bg-blue-50/50 px-3 py-2">
              <MousePointerClick className="h-4 w-4 text-blue-600" />
              <span className="text-sm font-medium text-blue-900">Actions</span>
            </div>
            <div className="flex-1 overflow-y-auto p-3">
              {mode === 'transfer' ? (
                <div className="space-y-3">
                  <button className="group flex w-full items-center justify-between rounded border border-slate-200 bg-white px-3 py-2 text-left text-xs font-medium text-slate-700 transition-colors hover:border-blue-400 hover:text-blue-600">
                    <span className="mr-2 truncate">Select Source Structure</span>
                    <span className="whitespace-nowrap rounded bg-slate-100 px-1.5 py-0.5 text-[10px] text-slate-500 group-hover:bg-blue-50 group-hover:text-blue-600">
                      None
                    </span>
                  </button>
                  <button className="group flex w-full items-center justify-between rounded border border-slate-200 bg-white px-3 py-2 text-left text-xs font-medium text-slate-700 transition-colors hover:border-blue-400 hover:text-blue-600">
                    <span className="mr-2 truncate">Select Target Structure</span>
                    <span className="whitespace-nowrap rounded bg-slate-100 px-1.5 py-0.5 text-[10px] text-slate-500 group-hover:bg-blue-50 group-hover:text-blue-600">
                      None
                    </span>
                  </button>
                  <div className="pt-2">
                    <button className="w-full rounded bg-emerald-600 py-2 text-sm font-medium text-white shadow-sm transition-colors hover:bg-emerald-700">
                      Apply Transfer
                    </button>
                  </div>
                </div>
              ) : null}

              {mode === 'supercell' ? (
                <div className="space-y-4">
                  <button className="w-full rounded border border-slate-200 bg-white px-3 py-2 text-xs font-medium text-slate-700 shadow-sm transition-colors hover:border-blue-400 hover:text-blue-600">
                    Select Base Structure
                  </button>

                  <div className="space-y-2 border-t border-slate-100 pt-2">
                    <div className="flex items-center justify-between text-xs text-slate-500">
                      <span>Supercell Grid</span>
                      <span className="rounded-full bg-slate-100 px-2 py-0.5 text-[10px]">
                        5x5
                      </span>
                    </div>
                    <div className="grid grid-cols-5 gap-1.5 rounded border border-slate-100 bg-slate-50 p-1">
                      {Array.from({ length: 25 }).map((_, i) => (
                        <button
                          key={i}
                          className={`aspect-square rounded-sm border transition-all duration-200 ${
                            i === 12
                              ? 'border-indigo-600 bg-indigo-500 shadow-sm'
                              : 'border-slate-200 bg-white hover:border-indigo-300'
                          }`}
                        />
                      ))}
                    </div>
                  </div>

                  <button className="mt-auto w-full rounded bg-indigo-600 py-2 text-sm font-medium text-white shadow-sm transition-colors hover:bg-indigo-700">
                    Build Supercell
                  </button>
                </div>
              ) : null}

              {mode === 'vibration' ? (
                <div className="flex h-full flex-col items-center justify-center gap-4">
                  <div className="rounded-full bg-rose-50 p-4">
                    <Activity className="h-8 w-8 text-rose-400" />
                  </div>
                  <div className="space-y-1 text-center">
                    <p className="text-sm font-medium text-slate-700">
                      Preview Only
                    </p>
                    <p className="rounded bg-slate-100 px-2 py-0.5 text-xs font-mono text-slate-400">
                      mode: phonopy
                    </p>
                  </div>
                  <button className="mt-2 flex w-full items-center justify-center gap-2 rounded-md bg-rose-600 py-2.5 font-medium text-white shadow-md transition-colors hover:bg-rose-700">
                    <Play className="h-4 w-4 fill-current" />
                    Run Preview
                  </button>
                </div>
              ) : null}
            </div>
          </div>
        </div>
      </div>
    </div>
  )
}
