import { Layers } from 'lucide-react'

import type { MolstarStructure } from './types'

import MolstarViewer from '@/components/molstar/MolstarViewer'
import { Button } from '@/components/ui/button'
import { Card, CardContent, CardHeader } from '@/components/ui/card'

export type ViewerPanelProps = {
  structures: Array<MolstarStructure>
}

export function ViewerPanel({ structures }: ViewerPanelProps) {
  return (
    <Card className="flex flex-1 flex-col border-white/10 bg-white/5 shadow-lg shadow-black/40">
      <CardHeader className="flex items-center justify-between">
        <div>
          <p className="text-xs uppercase tracking-[0.3em] text-white/50">
            3D Viewer
          </p>
          <p className="text-lg font-semibold">Mol* Preview</p>
        </div>
        <Layers className="h-4 w-4 text-white/40" />
      </CardHeader>
      <CardContent className="flex flex-1 flex-col gap-4">
        <div className="h-80 flex-1 lg:h-full">
          <MolstarViewer structures={structures} />
        </div>
        <div className="grid grid-cols-2 gap-2 text-xs text-white/70">
          <Button variant="outline" size="sm" disabled>
            Auto Fit
          </Button>
          <Button variant="outline" size="sm" disabled>
            Reset View
          </Button>
          <Button variant="outline" size="sm" disabled>
            Snapshot
          </Button>
          <Button variant="outline" size="sm" disabled>
            Clip
          </Button>
        </div>
      </CardContent>
    </Card>
  )
}
