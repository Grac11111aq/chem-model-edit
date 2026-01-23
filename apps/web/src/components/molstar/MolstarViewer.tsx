import { useEffect, useMemo, useRef, useState } from 'react'
import type { Viewer } from 'molstar/lib/apps/viewer/app'

import { cn } from '@/lib/utils'

type MolstarViewerStructure = {
  id: string
  pdbText?: string
  bcifUrl?: string
  opacity?: number
  visible?: boolean
}

type MolstarViewerProps = {
  pdbText?: string
  bcifUrl?: string
  structures?: Array<MolstarViewerStructure>
  onError?: (message: string) => void
  onLoad?: () => void
  selectedAtomIndices?: Array<number>
  disabledAtomIndices?: Array<number>
  onAtomToggle?: (index: number) => void
  className?: string
}

type MolstarHelpers = {
  StructureElement: any
  StructureProperties: any
  Bond: any
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
  selectedAtomIndices,
  disabledAtomIndices,
  onAtomToggle,
  className,
}: MolstarViewerProps) {
  const containerRef = useRef<HTMLDivElement | null>(null)
  const viewerRef = useRef<Viewer | null>(null)
  const lastSignatureRef = useRef<string | null>(null)
  const timerRef = useRef<number | null>(null)
  const activeRef = useRef(true)
  const abortRef = useRef<AbortController | null>(null)
  const structureReadyRef = useRef<string | null>(null)
  const selectionSignatureRef = useRef<string | null>(null)
  const pendingSelectionRef = useRef<Array<number> | null>(null)
  const latestSelectionRef = useRef<Array<number>>([])
  const helpersRef = useRef<MolstarHelpers | null>(null)
  const clickSubRef = useRef<{ unsubscribe: () => void } | null>(null)
  const onAtomToggleRef = useRef<typeof onAtomToggle>(onAtomToggle)
  const disabledSetRef = useRef<Set<number>>(new Set())
  const [viewerReady, setViewerReady] = useState(false)

  const normalizedStructures = useMemo<Array<MolstarViewerStructure>>(() => {
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

  const selectionEnabled = selectedAtomIndices !== undefined

  const normalizedSelection = useMemo(() => {
    if (!selectedAtomIndices || selectedAtomIndices.length === 0) {
      return []
    }
    const disabledSet = new Set(disabledAtomIndices ?? [])
    return Array.from(new Set(selectedAtomIndices))
      .filter((index) => !disabledSet.has(index))
      .sort((a, b) => a - b)
  }, [selectedAtomIndices, disabledAtomIndices])

  const selectionSignature = useMemo(
    () => normalizedSelection.join(','),
    [normalizedSelection],
  )

  useEffect(() => {
    if (!selectionEnabled) {
      return
    }
    latestSelectionRef.current = normalizedSelection
  }, [normalizedSelection, selectionEnabled])

  useEffect(() => {
    onAtomToggleRef.current = onAtomToggle
  }, [onAtomToggle])

  useEffect(() => {
    disabledSetRef.current = new Set(disabledAtomIndices ?? [])
  }, [disabledAtomIndices])

  const applySelection = (indices: Array<number>) => {
    const viewer = viewerRef.current as any
    if (!viewer || typeof viewer.structureInteractivity !== 'function') {
      return
    }
    const disabledSet = disabledSetRef.current
    const filtered =
      disabledSet.size > 0
        ? indices.filter((index) => !disabledSet.has(index))
        : indices
    try {
      viewer.structureInteractivity({ action: 'select' })
    } catch (error) {
      const message =
        error instanceof Error && error.message
          ? error.message
          : 'Mol* selection update failed.'
      console.error('Mol* selection clear failed.', error)
      onError?.(message)
      return
    }
    if (filtered.length === 0) {
      return
    }
    try {
      viewer.structureInteractivity({
        action: 'select',
        applyGranularity: false,
        elements: {
          items: { atom_index: filtered },
        },
      })
    } catch (error) {
      const message =
        error instanceof Error && error.message
          ? error.message
          : 'Mol* selection update failed.'
      console.error('Mol* selection apply failed.', error)
      onError?.(message)
    }
  }

  const loadBallAndStick = async (
    viewer: Viewer,
    items: Array<MolstarViewerStructure>,
    loadSignature: string,
  ) => {
    const plugin = (viewer as unknown as { plugin: any }).plugin
    abortRef.current?.abort()
    const controller = new AbortController()
    abortRef.current = controller
    structureReadyRef.current = null
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
          const structure =
            await plugin.builders.structure.createStructure(model)
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
        structureReadyRef.current = loadSignature
        if (selectionEnabled) {
          if (pendingSelectionRef.current) {
            const selection = pendingSelectionRef.current
            pendingSelectionRef.current = null
            applySelection(selection)
            selectionSignatureRef.current = selection.join(',')
          } else {
            applySelection(latestSelectionRef.current)
            selectionSignatureRef.current = latestSelectionRef.current.join(',')
          }
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
    abortRef.current?.abort()
    abortRef.current = null
    const plugin = (viewer as unknown as { plugin: any }).plugin
    await plugin.clear()
    structureReadyRef.current = null
    selectionSignatureRef.current = null
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
      clickSubRef.current?.unsubscribe()
      clickSubRef.current = null
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
      void loadBallAndStick(
        viewerRef.current,
        normalizedStructures,
        signature,
      ).then(() => {
        lastSignatureRef.current = signature
      })
    }, 150)
    return () => {
      if (timerRef.current) {
        window.clearTimeout(timerRef.current)
      }
    }
  }, [normalizedStructures, signature, viewerReady])

  useEffect(() => {
    if (lastSignatureRef.current !== signature) {
      selectionSignatureRef.current = null
      pendingSelectionRef.current = null
    }
  }, [signature])

  useEffect(() => {
    if (!viewerRef.current || !viewerReady) {
      return
    }
    if (selectionEnabled) {
      return
    }
    selectionSignatureRef.current = null
    pendingSelectionRef.current = null
    applySelection([])
  }, [selectionEnabled, viewerReady])

  useEffect(() => {
    if (!viewerRef.current || !viewerReady) {
      return
    }
    if (!selectionEnabled) {
      return
    }
    if (structureReadyRef.current !== signature) {
      pendingSelectionRef.current = normalizedSelection
      return
    }
    if (selectionSignatureRef.current === selectionSignature) {
      return
    }
    applySelection(normalizedSelection)
    selectionSignatureRef.current = selectionSignature
  }, [
    normalizedSelection,
    selectionSignature,
    selectionEnabled,
    signature,
    viewerReady,
  ])

  useEffect(() => {
    if (!viewerRef.current || !viewerReady) {
      return
    }
    if (!selectionEnabled) {
      clickSubRef.current?.unsubscribe()
      clickSubRef.current = null
      return
    }
    let cancelled = false
    const plugin = (viewerRef.current as unknown as { plugin: any }).plugin
    if (!plugin?.behaviors?.interaction?.click?.subscribe) {
      return
    }

    const setup = async () => {
      if (!helpersRef.current) {
        const [{ StructureElement }, { StructureProperties }, { Bond }] =
          await Promise.all([
            import('molstar/lib/mol-model/structure/structure/element'),
            import('molstar/lib/mol-model/structure/structure/properties'),
            import('molstar/lib/mol-model/structure/structure/unit/bonds'),
          ])
        helpersRef.current = { StructureElement, StructureProperties, Bond }
      }
      if (cancelled) {
        return
      }
      clickSubRef.current?.unsubscribe()
      clickSubRef.current = plugin.behaviors.interaction.click.subscribe(
        (event: any) => {
          const toggle = onAtomToggleRef.current
          if (!toggle) {
            return
          }
          const helpers = helpersRef.current
          if (!helpers) {
            return
          }
          const loci = event?.current
          let elementLoci = loci
          if (helpers.Bond.isLoci(loci)) {
            elementLoci = helpers.Bond.toFirstStructureElementLoci(loci)
          }
          if (!helpers.StructureElement.Loci.is(elementLoci)) {
            return
          }
          const location =
            helpers.StructureElement.Loci.getFirstLocation(elementLoci)
          if (!location) {
            return
          }
          const index = helpers.StructureProperties.atom.sourceIndex(location)
          if (typeof index !== 'number' || Number.isNaN(index)) {
            return
          }
          if (disabledSetRef.current.has(index)) {
            window.setTimeout(() => {
              applySelection(latestSelectionRef.current)
            }, 0)
            return
          }
          toggle(index)
        },
      )
    }

    void setup()

    return () => {
      cancelled = true
      clickSubRef.current?.unsubscribe()
      clickSubRef.current = null
    }
  }, [selectionEnabled, viewerReady])

  return (
    <div
      ref={containerRef}
      className={cn(
        'relative h-full w-full overflow-hidden rounded-xl border border-white/10',
        className,
      )}
    />
  )
}
