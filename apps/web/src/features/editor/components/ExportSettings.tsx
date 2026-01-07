import { Badge } from '@/components/ui/badge'
import { Card, CardContent, CardHeader } from '@/components/ui/card'

export function ExportSettings() {
  return (
    <Card className="border-white/10 bg-white/5">
      <CardHeader className="flex items-center justify-between">
        <p className="text-xs uppercase tracking-[0.3em] text-white/50">
          Export Settings
        </p>
        <span className="text-xs text-white/40">Fixed (v1)</span>
      </CardHeader>
      <CardContent className="space-y-3 text-sm text-white/70">
        <div className="flex items-center justify-between">
          <span>Format</span>
          <Badge
            variant="outline"
            className="border-white/10 px-3 py-1 text-xs"
          >
            Match Input
          </Badge>
        </div>
        <div className="flex items-center justify-between">
          <span>Units</span>
          <Badge
            variant="outline"
            className="border-white/10 px-3 py-1 text-xs"
          >
            Match Input
          </Badge>
        </div>
        <div className="flex items-center justify-between">
          <span>Coordinates</span>
          <Badge
            variant="outline"
            className="border-white/10 px-3 py-1 text-xs"
          >
            Match Input
          </Badge>
        </div>
        <div className="flex items-center justify-between">
          <span>celldm</span>
          <Badge
            variant="outline"
            className="border-white/10 px-3 py-1 text-xs"
          >
            Keep
          </Badge>
        </div>
        <p className="text-xs text-white/40">
          出力は入力形式・単位を保持します。将来のオプション拡張に備えた枠です。
        </p>
      </CardContent>
    </Card>
  )
}
