import { RefreshCw } from 'lucide-react'

import type { Atom } from '@/lib/types'
import type { PbcState, StructureState } from '../types'

import { Badge } from '@/components/ui/badge'
import { Button } from '@/components/ui/button'
import { Card, CardContent, CardHeader } from '@/components/ui/card'
import {
  Select,
  SelectContent,
  SelectItem,
  SelectTrigger,
  SelectValue,
} from '@/components/ui/select'
import { cn } from '@/lib/utils'

type DistanceRow = {
  index: number
  source: Atom
  target: Atom
  dist: number
}

type ComparePanelProps = {
  compareCandidates: Array<StructureState>
  compareTarget: StructureState | null
  countMismatch: boolean
  pbcState: PbcState
  distanceRows: Array<DistanceRow>
  selectedIndices: Array<number>
  atomsCount: number
  selectedDistance: number | null
  activeId: string
  onChangeTarget: (id: string | null) => void
  onAlignSelected: () => void
  onSurfaceTransfer: () => void
}

export function ComparePanel({
  compareCandidates,
  compareTarget,
  countMismatch,
  pbcState,
  distanceRows,
  selectedIndices,
  atomsCount,
  selectedDistance,
  activeId,
  onChangeTarget,
  onAlignSelected,
  onSurfaceTransfer,
}: ComparePanelProps) {
  return (
    <Card className="border-white/10 bg-white/5">
      <CardHeader className="space-y-3">
        <div className="flex items-center justify-between">
          <p className="text-xs uppercase tracking-[0.3em] text-white/50">
            Compare
          </p>
          <RefreshCw className="h-4 w-4 text-white/40" />
        </div>
      </CardHeader>
      <CardContent className="space-y-4 text-sm text-white/70">
        <div className="flex items-center justify-between gap-3">
          <span>Compare With</span>
          <Select
            value={compareTarget?.id ?? ''}
            onValueChange={(value) => onChangeTarget(value || null)}
          >
            <SelectTrigger className="w-[160px] border-white/10 bg-slate-950/70 text-xs text-white/80">
              <SelectValue placeholder="None" />
            </SelectTrigger>
            <SelectContent>
              {compareCandidates.length === 0 ? (
                <SelectItem value="">None</SelectItem>
              ) : null}
              {compareCandidates.map((candidate) => (
                <SelectItem key={candidate.id} value={candidate.id}>
                  {candidate.name}
                </SelectItem>
              ))}
            </SelectContent>
          </Select>
        </div>

        {countMismatch ? (
          <div className="rounded-lg border border-rose-400/30 bg-rose-500/10 px-3 py-2 text-xs text-rose-200">
            原子数が一致しません。インデックス対応の転写・距離に影響します。
          </div>
        ) : null}

        <div className="flex items-center justify-between">
          <span>PBC</span>
          <Badge
            variant="outline"
            className={cn(
              'rounded-full border-white/10 px-3 py-1 text-xs',
              pbcState.enabled ? 'text-emerald-200' : 'text-white/60',
            )}
          >
            {pbcState.enabled ? 'On' : 'Off'}
          </Badge>
        </div>
        {!pbcState.enabled && compareTarget ? (
          <div className="rounded-lg border border-white/10 bg-slate-900/60 px-3 py-2 text-xs text-white/50">
            {pbcState.reason}
          </div>
        ) : null}

        <div className="flex items-center justify-between">
          <span>RMSD</span>
          <span className="font-semibold text-white">0.042 Å</span>
        </div>
        <div className="flex items-center justify-between">
          <span>Selected Pair</span>
          <span className="font-semibold text-white">
            {selectedDistance ? `${selectedDistance.toFixed(4)} Å` : '—'}
          </span>
        </div>

        <div className="rounded-xl border border-white/10 bg-slate-900/60 px-3 py-2 text-xs text-white/60">
          <p className="uppercase tracking-[0.3em] text-white/40">
            A↔B Distances (Index)
          </p>
          <div className="mt-2 space-y-1">
            {compareTarget ? null : (
              <p className="text-white/40">比較対象を選択してください。</p>
            )}
            {compareTarget && distanceRows.length === 0 ? (
              <p className="text-white/40">距離を表示できません。</p>
            ) : null}
            {compareTarget
              ? distanceRows.map((row) => (
                  <div
                    key={`distance-${row.index}`}
                    className="flex justify-between"
                  >
                    <span>
                      #{row.index + 1} {row.source.symbol}→{row.target.symbol}
                    </span>
                    <span className="text-white">{row.dist.toFixed(4)} Å</span>
                  </div>
                ))
              : null}
          </div>
        </div>

        <div className="flex items-center justify-between">
          <span>Selected Atoms</span>
          <span className="font-semibold text-white">
            {selectedIndices.length} / {atomsCount}
          </span>
        </div>

        <Button
          variant="secondary"
          size="sm"
          className="w-full bg-white/10 text-xs text-white hover:bg-white/20"
          onClick={onAlignSelected}
        >
          Align Selected (Centroid)
        </Button>
        <Button
          variant="outline"
          size="sm"
          className="w-full border-amber-200/40 text-xs text-amber-100/90 hover:text-amber-200"
          onClick={onSurfaceTransfer}
          disabled={!compareTarget}
        >
          Surface Transfer ({activeId}→{compareTarget?.id ?? '—'})
        </Button>
      </CardContent>
    </Card>
  )
}
