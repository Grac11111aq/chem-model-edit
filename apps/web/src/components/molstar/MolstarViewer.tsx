import { useEffect, useRef } from 'react'
import { Viewer } from 'molstar/lib/apps/viewer/app'
import 'molstar/lib/mol-plugin-ui/skin/light.scss'

type MolstarViewerProps = {
  pdbText?: string
}

export default function MolstarViewer({ pdbText }: MolstarViewerProps) {
  const containerRef = useRef<HTMLDivElement | null>(null)
  const viewerRef = useRef<Viewer | null>(null)
  const lastPdbRef = useRef<string | null>(null)
  const timerRef = useRef<number | null>(null)

  const loadBallAndStick = async (viewer: Viewer, pdb: string) => {
    const plugin = (viewer as unknown as { plugin: any }).plugin
    await plugin.clear()
    const data = await plugin.builders.data.rawData({ data: pdb }, { state: { isGhost: true } })
    const trajectory = await plugin.builders.structure.parseTrajectory(data, 'pdb')
    const model = await plugin.builders.structure.createModel(trajectory)
    const structure = await plugin.builders.structure.createStructure(model)
    const component = await plugin.builders.structure.tryCreateComponentStatic(structure, 'all')
    if (component) {
      await plugin.builders.structure.representation.addRepresentation(component, {
        type: 'ball-and-stick',
        color: 'element-symbol',
      })
    }
    plugin.managers?.camera?.reset?.()
  }

  useEffect(() => {
    if (!containerRef.current || viewerRef.current) {
      return
    }

    let active = true

    Viewer.create(containerRef.current, {
      layoutIsExpanded: false,
      layoutShowControls: false,
      layoutShowLeftPanel: false,
      layoutShowSequence: false,
      layoutShowLog: false,
      layoutShowSidebar: false,
      viewportShowExpand: false,
      viewportShowControls: false,
      backgroundColor: 0x0b1120,
    }).then((viewer) => {
      if (!active) {
        viewer.dispose()
        return
      }
      viewerRef.current = viewer
    })

    return () => {
      active = false
      if (viewerRef.current) {
        viewerRef.current.dispose()
        viewerRef.current = null
      }
    }
  }, [])

  useEffect(() => {
    if (!viewerRef.current || !pdbText) {
      return
    }
    if (lastPdbRef.current === pdbText) {
      return
    }
    if (timerRef.current) {
      window.clearTimeout(timerRef.current)
    }
    timerRef.current = window.setTimeout(() => {
      if (!viewerRef.current) {
        return
      }
      void loadBallAndStick(viewerRef.current, pdbText).then(() => {
        lastPdbRef.current = pdbText
      })
    }, 150)
    return () => {
      if (timerRef.current) {
        window.clearTimeout(timerRef.current)
      }
    }
  }, [pdbText])

  return (
    <div
      ref={containerRef}
      className="h-full w-full overflow-hidden rounded-xl border border-white/10"
    />
  )
}
