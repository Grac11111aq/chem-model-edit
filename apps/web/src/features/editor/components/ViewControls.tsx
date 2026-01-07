import { Droplets } from 'lucide-react'

import { Badge } from '@/components/ui/badge'
import { Card, CardContent, CardHeader } from '@/components/ui/card'

type ViewControlsProps = {
  overlayEnabled: boolean
  visibleOverlayCount: number
  activeOpacity: number
}

export function ViewControls({
  overlayEnabled,
  visibleOverlayCount,
  activeOpacity,
}: ViewControlsProps) {
  return (
    <Card className="border-white/10 bg-white/5">
      <CardHeader className="flex items-center justify-between">
        <p className="text-xs uppercase tracking-[0.3em] text-white/50">
          View Controls
        </p>
        <Droplets className="h-4 w-4 text-white/40" />
      </CardHeader>
      <CardContent className="space-y-3 text-sm text-white/70">
        <div className="flex items-center justify-between">
          <span>Overlay Mode</span>
          <Badge
            variant="outline"
            className="border-white/10 px-3 py-1 text-xs"
          >
            {overlayEnabled ? 'On' : 'Off'}
          </Badge>
        </div>
        <div className="flex items-center justify-between">
          <span>Visible Structures</span>
          <Badge
            variant="outline"
            className="border-white/10 px-3 py-1 text-xs"
          >
            {visibleOverlayCount}
          </Badge>
        </div>
        <div className="flex items-center justify-between">
          <span>Active Opacity</span>
          <span className="text-xs">{Math.round(activeOpacity * 100)}%</span>
        </div>
      </CardContent>
    </Card>
  )
}
