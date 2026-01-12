import { useEffect, useMemo, useRef, useState } from 'react'
import type { Viewer } from 'molstar/lib/apps/viewer/app'

type MolstarStructure = {
  id: string
  pdbText?: string
  bcifUrl?: string
  opacity?: number
  visible?: boolean
}

type MolstarViewerProps = {
  pdbText?: string
  bcifUrl?: string
  structures?: Array<MolstarStructure>
  onError?: (message: string) => void
  onLoad?: () => void
}

const hashString = (value: string) => {
  let hash = 0
  for (let i = 0; i < value.length; i += 1) {
    hash = (hash << 5) - hash + value.charCodeAt(i)
    hash |= 0
  }
  return hash.toString(16)
}

const signatureForText = (value: string) => {
  const head = value.slice(0, 64)
  const tail = value.slice(-64)
  return `${value.length}:${hashString(value)}:${head}:${tail}`
}

export default function MolstarViewer({
  pdbText,
  bcifUrl,
  structures,
  onError,
  onLoad,
}: MolstarViewerProps) {
  const containerRef = useRef<HTMLDivElement | null>(null)
  const viewerRef = useRef<Viewer | null>(null)
  const lastSignatureRef = useRef<string | null>(null)
  const timerRef = useRef<number | null>(null)
  const activeRef = useRef(true)
  const abortRef = useRef<AbortController | null>(null)
  const [viewerReady, setViewerReady] = useState(false)

  const normalizedStructures = useMemo<Array<MolstarStructure>>(() => {
    if (structures && structures.length > 0) {
      return structures
    }
    if (bcifUrl) {
      return [{ id: 'single', bcifUrl, opacity: 1, visible: true }]
    }
    if (pdbText) {
      return [{ id: 'single', pdbText, opacity: 1, visible: true }]
    }
    return []
  }, [bcifUrl, pdbText, structures])

  const signature = useMemo(() => {
    if (normalizedStructures.length === 0) {
      return ''
    }
    return normalizedStructures
      .map((structure) => {
        const opacity = structure.opacity ?? 1
        const visible = structure.visible ?? true
        const payloadSignature = structure.bcifUrl
          ? `bcif:${structure.bcifUrl}`
          : structure.pdbText
            ? `pdb:${signatureForText(structure.pdbText)}`
            : 'empty'
        return `${structure.id}|${opacity}|${visible}|${payloadSignature}`
      })
      .join('::')
  }, [normalizedStructures])

  const loadBallAndStick = async (
    viewer: Viewer,
    items: Array<MolstarStructure>,
  ) => {
    const plugin = (viewer as unknown as { plugin: any }).plugin
    abortRef.current?.abort()
    const controller = new AbortController()
    abortRef.current = controller
    const isCancelled = () =>
      !activeRef.current || abortRef.current !== controller
    try {
      if (isCancelled()) {
        return
      }
      await plugin.clear()
      let loadedCount = 0
      let lastError: string | null = null
      for (const item of items) {
        if (isCancelled()) {
          break
        }
        if (item.visible === false) {
          continue
        }
        try {
          let data = null
          let format: 'pdb' | 'bcif' = 'pdb'
          if (item.bcifUrl) {
            const response = await fetch(item.bcifUrl, {
              signal: controller.signal,
            })
            if (!response.ok) {
              lastError = `BCIF fetch failed (HTTP ${response.status})`
              console.error('Mol* bcif 取得に失敗しました。', response.status)
              continue
            }
            if (isCancelled()) {
              break
            }
            const buffer = await response.arrayBuffer()
            if (isCancelled()) {
              break
            }
            data = await plugin.builders.data.rawData(
              { data: new Uint8Array(buffer) },
              { state: { isGhost: true } },
            )
            format = 'bcif'
          } else if (item.pdbText) {
            data = await plugin.builders.data.rawData(
              { data: item.pdbText },
              { state: { isGhost: true } },
            )
            format = 'pdb'
          }

          if (!data) {
            lastError = 'No structure data was provided.'
            continue
          }

          const trajectory = await plugin.builders.structure.parseTrajectory(
            data,
            format,
          )
          if (isCancelled()) {
            break
          }
          const model = await plugin.builders.structure.createModel(trajectory)
          if (isCancelled()) {
            break
          }
          const structure = await plugin.builders.structure.createStructure(
            model,
          )
          if (isCancelled()) {
            break
          }
          const component =
            await plugin.builders.structure.tryCreateComponentStatic(
              structure,
              'all',
            )
          if (component) {
            const opacity = Math.min(1, Math.max(0, item.opacity ?? 1))
            await plugin.builders.structure.representation.addRepresentation(
              component,
              {
                type: 'ball-and-stick',
                color: 'element-symbol',
                typeParams: { alpha: opacity },
              },
            )
            loadedCount += 1
          }
        } catch (error) {
          if (isCancelled()) {
            break
          }
          const message =
            error instanceof Error && error.message
              ? error.message
              : 'Mol* failed to render structure.'
          lastError = message
          console.error('Mol* Viewer load failed.', error)
        }
      }
      if (activeRef.current && abortRef.current === controller) {
        if (loadedCount > 0) {
          onLoad?.()
        } else if (lastError) {
          onError?.(lastError)
        }
        plugin.managers?.camera?.reset?.()
      }
    } finally {
      if (abortRef.current === controller) {
        abortRef.current = null
      }
    }
  }

  const clearViewer = async (viewer: Viewer) => {
    const plugin = (viewer as unknown as { plugin: any }).plugin
    await plugin.clear()
  }

  useEffect(() => {
    const container = containerRef.current
    if (!container || viewerRef.current) {
      return
    }

    activeRef.current = true

    void (async () => {
      try {
        await import('molstar/lib/mol-plugin-ui/skin/light.scss')
        const { Viewer } = await import('molstar/lib/apps/viewer/app')
        const { CifCoreProvider } =
          await import('molstar/lib/mol-plugin-state/formats/trajectory')
        const viewer = await Viewer.create(container, {
          layoutIsExpanded: false,
          layoutShowControls: false,
          layoutShowLeftPanel: false,
          layoutShowSequence: false,
          layoutShowLog: false,
          layoutShowSidebar: false,
          viewportShowExpand: false,
          viewportShowControls: false,
          backgroundColor: 0x0b1120,
          customFormats: [['bcif', CifCoreProvider]],
        })
        if (!activeRef.current) {
          viewer.dispose()
          return
        }
        viewerRef.current = viewer
        setViewerReady(true)
      } catch (error) {
        console.error('Mol* Viewerの初期化に失敗しました。', error)
      }
    })()

    return () => {
      activeRef.current = false
      abortRef.current?.abort()
      abortRef.current = null
      if (viewerRef.current) {
        viewerRef.current.dispose()
        viewerRef.current = null
      }
      setViewerReady(false)
    }
  }, [])

  useEffect(() => {
    if (!viewerRef.current || !viewerReady) {
      return
    }
    if (signature.length === 0) {
      if (lastSignatureRef.current !== signature) {
        void clearViewer(viewerRef.current).then(() => {
          lastSignatureRef.current = signature
        })
      }
      return
    }
    if (lastSignatureRef.current === signature) {
      return
    }
    if (timerRef.current) {
      window.clearTimeout(timerRef.current)
    }
    timerRef.current = window.setTimeout(() => {
      if (!viewerRef.current) {
        return
      }
      void loadBallAndStick(viewerRef.current, normalizedStructures).then(
        () => {
          lastSignatureRef.current = signature
        },
      )
    }, 150)
    return () => {
      if (timerRef.current) {
        window.clearTimeout(timerRef.current)
      }
    }
  }, [normalizedStructures, signature, viewerReady])

  return (
    <div
      ref={containerRef}
      className="relative h-full w-full overflow-hidden rounded-xl border border-white/10"
    />
  )
}
