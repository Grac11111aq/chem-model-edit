import { PencilRuler } from 'lucide-react'

import { Button } from '@/components/ui/button'
import { Card, CardContent, CardHeader } from '@/components/ui/card'
import { Input } from '@/components/ui/input'

const axes: Array<'x' | 'y' | 'z'> = ['x', 'y', 'z']

type SelectionToolsProps = {
  shiftDraft: { x: string; y: string; z: string }
  onShiftDraftChange: (axis: 'x' | 'y' | 'z', value: string) => void
  onApplyShift: () => void
  onAlignOrigin: () => void
  onAlignCentroid: () => void
  onCopySelected: () => void
  onPasteAppend: () => void
  onClearSelection: () => void
}

export function SelectionTools({
  shiftDraft,
  onShiftDraftChange,
  onApplyShift,
  onAlignOrigin,
  onAlignCentroid,
  onCopySelected,
  onPasteAppend,
  onClearSelection,
}: SelectionToolsProps) {
  return (
    <Card className="border-white/10 bg-white/5">
      <CardHeader className="flex items-center justify-between">
        <p className="text-xs uppercase tracking-[0.3em] text-white/50">
          Selection Tools
        </p>
        <PencilRuler className="h-4 w-4 text-white/40" />
      </CardHeader>
      <CardContent className="space-y-4">
        <div className="rounded-xl border border-white/10 bg-slate-900/60 p-3">
          <p className="text-white/60">Shift by vector</p>
          <div className="mt-2 grid grid-cols-3 gap-2 text-xs">
            {axes.map((axis) => (
              <Input
                key={axis}
                value={shiftDraft[axis]}
                onChange={(event) =>
                  onShiftDraftChange(axis, event.target.value)
                }
                className="h-9 border-white/10 bg-white/5 text-white/70 focus-visible:ring-amber-300"
                placeholder={`d${axis}`}
              />
            ))}
          </div>
        </div>
        <div className="grid grid-cols-2 gap-2 text-xs">
          <Button
            variant="secondary"
            size="sm"
            className="bg-white/10 text-white hover:bg-white/20"
            onClick={onApplyShift}
          >
            Apply Shift
          </Button>
          <Button
            variant="outline"
            size="sm"
            className="border-white/10 text-white/70 hover:text-white"
            onClick={onAlignOrigin}
          >
            Align to Origin
          </Button>
          <Button
            variant="outline"
            size="sm"
            className="border-white/10 text-white/70 hover:text-white"
            onClick={onAlignCentroid}
          >
            Align to Centroid
          </Button>
          <Button
            variant="outline"
            size="sm"
            className="border-white/10 text-white/70 hover:text-white"
            onClick={onCopySelected}
          >
            Copy Selected
          </Button>
          <Button
            variant="outline"
            size="sm"
            className="border-white/10 text-white/70 hover:text-white"
            onClick={onPasteAppend}
          >
            Paste Append
          </Button>
          <Button
            variant="outline"
            size="sm"
            className="border-white/10 text-white/70 hover:text-white"
            onClick={onClearSelection}
          >
            Clear Selection
          </Button>
        </div>
      </CardContent>
    </Card>
  )
}
