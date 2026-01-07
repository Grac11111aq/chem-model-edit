import { Fragment } from 'react'

import type { Atom } from '@/lib/types'
import type { AtomDraft } from '../types'

import { Button } from '@/components/ui/button'
import { Card, CardContent, CardHeader } from '@/components/ui/card'
import { Checkbox } from '@/components/ui/checkbox'
import { Input } from '@/components/ui/input'
import { Textarea } from '@/components/ui/textarea'
import { cn } from '@/lib/utils'

type AtomEditorPanelProps = {
  qeInput: string
  isParsing: boolean
  error: string | null
  exported: string
  atoms: Array<Atom>
  atomDrafts: Array<AtomDraft>
  selectedIndices: Array<number>
  onChangeQeInput: (value: string) => void
  onParse: () => void
  onResetSample: () => void
  onExport: () => void
  onShare: () => void
  onAddAtom: () => void
  onToggleSelect: (index: number) => void
  onAtomChange: (index: number, field: keyof AtomDraft, value: string) => void
}

export function AtomEditorPanel({
  qeInput,
  isParsing,
  error,
  exported,
  atoms,
  atomDrafts,
  selectedIndices,
  onChangeQeInput,
  onParse,
  onResetSample,
  onExport,
  onShare,
  onAddAtom,
  onToggleSelect,
  onAtomChange,
}: AtomEditorPanelProps) {
  return (
    <Card className="border-white/10 bg-white/5 shadow-lg shadow-black/40">
      <CardHeader className="flex flex-wrap items-center justify-between gap-3">
        <div>
          <p className="text-sm uppercase tracking-[0.3em] text-white/50">
            Atom Table
          </p>
          <p className="text-lg font-semibold">Editable Coordinates</p>
        </div>
        <div className="flex flex-wrap items-center gap-2 text-xs text-white/60">
          <span className="rounded-full border border-white/10 px-3 py-1">
            XYZ mode
          </span>
          <span className="rounded-full border border-white/10 px-3 py-1">
            Ångström
          </span>
          <Button
            variant="outline"
            size="sm"
            className="rounded-full border-white/10 text-white/70 hover:text-white"
            onClick={onExport}
          >
            Export .in
          </Button>
          <Button
            variant="outline"
            size="sm"
            className="rounded-full border-amber-200/40 text-amber-100/90 hover:text-amber-200"
            onClick={onShare}
          >
            Share HTML
          </Button>
        </div>
      </CardHeader>
      <CardContent className="space-y-4">
        <div className="rounded-xl border border-white/10 bg-slate-900/60 p-4">
          <div className="flex flex-wrap items-center justify-between gap-2 text-xs text-white/60">
            <span className="uppercase tracking-[0.3em]">Input</span>
            <div className="flex items-center gap-2">
              <Button
                variant="outline"
                size="sm"
                className="rounded-full border-white/10 text-white/70 hover:text-white"
                onClick={onResetSample}
              >
                Reset
              </Button>
              <Button
                size="sm"
                className="rounded-full bg-amber-300 text-slate-900 hover:bg-amber-200"
                onClick={onParse}
                disabled={isParsing}
              >
                {isParsing ? 'Parsing...' : 'Parse'}
              </Button>
            </div>
          </div>
          <Textarea
            value={qeInput}
            onChange={(event) => onChangeQeInput(event.target.value)}
            className="mt-3 h-32 resize-none border-white/10 bg-slate-950/80 font-mono text-xs text-white/80 focus-visible:ring-amber-300"
            aria-label="Quantum ESPRESSO input"
          />
          {error ? <p className="mt-2 text-xs text-rose-300">{error}</p> : null}
          {exported ? (
            <p className="mt-2 text-xs text-emerald-200">
              エクスポート済み（クリップボードにコピーしました）
            </p>
          ) : null}
        </div>

        <div className="overflow-hidden rounded-xl border border-white/10">
          <div className="flex items-center justify-between border-b border-white/10 bg-slate-900/70 px-3 py-2 text-xs text-white/60">
            <span>Atoms</span>
            <Button
              variant="outline"
              size="sm"
              className="rounded-full border-white/10 text-xs text-white/70 hover:text-white"
              onClick={onAddAtom}
            >
              Add Atom
            </Button>
          </div>
          <div className="max-h-[360px] overflow-auto bg-slate-950/60">
            <table className="w-full text-left text-xs">
              <thead className="sticky top-0 bg-slate-950/90 text-white/60">
                <tr>
                  {['Sel', '#', 'Atom', 'x', 'y', 'z'].map((label) => (
                    <th key={label} className="px-3 py-2 font-medium">
                      {label}
                    </th>
                  ))}
                </tr>
              </thead>
              <tbody>
                {atoms.length === 0 ? (
                  <tr>
                    <td
                      colSpan={6}
                      className="px-3 py-6 text-center text-white/50"
                    >
                      まだ構造が読み込まれていません
                    </td>
                  </tr>
                ) : null}
                {atoms.map((atom, index) => {
                  const draft =
                    atomDrafts[index] ??
                    ({
                      symbol: atom.symbol,
                      x: atom.x.toFixed(4),
                      y: atom.y.toFixed(4),
                      z: atom.z.toFixed(4),
                    } satisfies AtomDraft)
                  return (
                    <Fragment key={`${atom.symbol}-${index}`}>
                      <tr className="border-t border-white/5">
                        <td className="px-3 py-2">
                          <Checkbox
                            checked={selectedIndices.includes(index)}
                            onCheckedChange={() => onToggleSelect(index)}
                            className={cn(
                              'h-4 w-4 rounded border-white/20 data-[state=checked]:bg-amber-300 data-[state=checked]:text-slate-900',
                            )}
                            aria-label={`select atom ${index + 1}`}
                          />
                        </td>
                        <td className="px-3 py-2 text-white/40">{index + 1}</td>
                        <td className="px-3 py-2">
                          <Input
                            value={draft.symbol}
                            onChange={(event) =>
                              onAtomChange(index, 'symbol', event.target.value)
                            }
                            className="h-8 w-16 border-white/10 bg-slate-900/70 text-white/80 focus-visible:ring-amber-300"
                            aria-label={`atom ${index + 1} symbol`}
                          />
                        </td>
                        {(['x', 'y', 'z'] as const).map((axis) => (
                          <td key={axis} className="px-3 py-2">
                            <Input
                              value={draft[axis]}
                              onChange={(event) =>
                                onAtomChange(index, axis, event.target.value)
                              }
                              className="h-8 w-24 border-white/10 bg-slate-900/70 text-white/80 focus-visible:ring-amber-300"
                              aria-label={`atom ${index + 1} ${axis}`}
                            />
                          </td>
                        ))}
                      </tr>
                    </Fragment>
                  )
                })}
              </tbody>
            </table>
          </div>
        </div>
      </CardContent>
    </Card>
  )
}
