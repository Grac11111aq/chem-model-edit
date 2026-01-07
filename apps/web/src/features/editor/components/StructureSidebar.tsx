import { ArrowDownToLine, ChevronDown, Plus } from 'lucide-react'
import type { RefObject } from 'react'

import type { StructureState } from '../types'

import { Button } from '@/components/ui/button'
import { Card, CardContent, CardHeader } from '@/components/ui/card'
import { Checkbox } from '@/components/ui/checkbox'
import { Slider } from '@/components/ui/slider'
import { cn } from '@/lib/utils'

const overlayButtonStyles = {
  on: 'border-primary/60 text-primary',
  off: 'border-border text-muted-foreground',
}

type StructureSidebarProps = {
  structures: Array<StructureState>
  activeId: string
  overlayEnabled: boolean
  fileInputRef: RefObject<HTMLInputElement | null>
  onToggleOverlay: () => void
  onAddStructure: () => void
  onSelectStructure: (id: string) => void
  onToggleVisibility: (id: string, next: boolean) => void
  onOpacityChange: (id: string, value: number) => void
  onImportFile: (event: React.ChangeEvent<HTMLInputElement>) => void
  onImportClipboard: () => void
}

export function StructureSidebar({
  structures,
  activeId,
  overlayEnabled,
  fileInputRef,
  onToggleOverlay,
  onAddStructure,
  onSelectStructure,
  onToggleVisibility,
  onOpacityChange,
  onImportFile,
  onImportClipboard,
}: StructureSidebarProps) {
  return (
    <Card className="border-white/10 bg-white/5 shadow-lg shadow-black/30">
      <CardHeader className="space-y-3">
        <div className="flex items-center justify-between">
          <p className="text-xs uppercase tracking-[0.3em] text-white/50">
            Structures
          </p>
          <div className="flex items-center gap-2">
            <Button
              variant="outline"
              size="sm"
              onClick={onToggleOverlay}
              className={cn(
                'rounded-full border px-3 text-[10px] uppercase tracking-[0.3em]',
                overlayEnabled
                  ? overlayButtonStyles.on
                  : overlayButtonStyles.off,
              )}
            >
              Overlay {overlayEnabled ? 'On' : 'Off'}
            </Button>
            <Button
              variant="outline"
              size="icon"
              onClick={onAddStructure}
              className="h-8 w-8 rounded-full border-white/10 bg-white/5 text-white/70 hover:text-white"
            >
              <Plus className="h-4 w-4" />
            </Button>
          </div>
        </div>
      </CardHeader>
      <CardContent className="space-y-4">
        <div className="space-y-3">
          {structures.map((structure) => (
            <div
              key={structure.id}
              className={cn(
                'group cursor-pointer rounded-xl border bg-slate-900/60 p-3 transition',
                structure.id === activeId
                  ? 'border-amber-200/60 shadow-lg shadow-amber-400/20'
                  : 'border-white/10 hover:border-white/30',
              )}
              onClick={() => onSelectStructure(structure.id)}
            >
              <div className="flex items-center justify-between">
                <div className="flex items-center gap-3">
                  <div
                    className={cn(
                      'h-10 w-10 rounded-xl bg-gradient-to-br',
                      structure.color,
                    )}
                  />
                  <div>
                    <p className="text-sm font-semibold">{structure.name}</p>
                    <p className="text-xs text-white/50">ID {structure.id}</p>
                  </div>
                </div>
                <ChevronDown className="h-4 w-4 text-white/40 transition group-hover:text-white" />
              </div>
              <div className="mt-3 flex items-center gap-2 text-xs text-white/60">
                <span className="rounded-full border border-white/10 px-2 py-1">
                  {structure.atoms.length} atoms
                </span>
                <span className="rounded-full border border-white/10 px-2 py-1">
                  {structure.id === activeId ? 'Active' : 'Idle'}
                </span>
              </div>
              <div
                className="mt-3 space-y-3 text-xs text-white/60"
                onClick={(event) => event.stopPropagation()}
              >
                <label className="flex items-center justify-between gap-2">
                  <span>Visible</span>
                  <Checkbox
                    checked={structure.isVisible}
                    onCheckedChange={(checked) =>
                      onToggleVisibility(structure.id, Boolean(checked))
                    }
                    className="border-white/20 data-[state=checked]:bg-amber-300 data-[state=checked]:text-slate-900"
                  />
                </label>
                <div className="flex items-center justify-between gap-3">
                  <span>Opacity</span>
                  <div className="flex w-full items-center gap-3">
                    <Slider
                      value={[Math.round(structure.opacity * 100)]}
                      min={0}
                      max={100}
                      step={1}
                      onValueChange={(value) =>
                        onOpacityChange(structure.id, value[0] / 100)
                      }
                      className="w-full"
                    />
                    <span className="w-10 text-right text-white/70">
                      {Math.round(structure.opacity * 100)}%
                    </span>
                  </div>
                </div>
              </div>
            </div>
          ))}
        </div>

        <div className="grid gap-2">
          <Button
            variant="outline"
            size="sm"
            className="w-full border-white/10 bg-white/5 text-xs text-white/70 hover:border-white/30"
            onClick={() => fileInputRef.current?.click()}
          >
            <ArrowDownToLine className="h-4 w-4" />
            Import .in File
          </Button>
          <Button
            variant="outline"
            size="sm"
            className="w-full border-white/10 bg-white/5 text-xs text-white/70 hover:border-white/30"
            onClick={onImportClipboard}
          >
            Import from Clipboard
          </Button>
          <input
            ref={fileInputRef}
            type="file"
            accept=".in,.txt"
            className="hidden"
            onChange={onImportFile}
          />
        </div>
      </CardContent>
    </Card>
  )
}
