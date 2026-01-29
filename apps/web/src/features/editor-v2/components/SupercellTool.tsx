import { useCallback, useEffect, useMemo, useRef, useState } from 'react'
import {
  AlertTriangle,
  Check,
  ChevronDown,
  Download,
  MoveDiagonal2,
  Sparkles,
} from 'lucide-react'

import type { CSSProperties, PointerEvent as ReactPointerEvent } from 'react'
import type { SupercellBuildMeta, SupercellGridAxis } from '@/lib/types'
import type { WorkspaceFile } from '../types'

import MolstarViewer from '@/components/molstar/MolstarViewer'
import { buildSupercell, getStructure, structureViewUrl } from '@/lib/api'
import { cn } from '@/lib/utils'
import { atomsToXyz } from '@/lib/xyz'

type SupercellResult = {
  structureId: string
  meta: SupercellBuildMeta
}

interface SupercellToolProps {
  structures: Array<WorkspaceFile>
  onSupercellCreated?: (result: SupercellResult) => void
}

type AxisMode = 'row-b' | 'row-a'

const DEFAULT_GRID_SIZE = 2
const GRID_CELL_PX = 42
const GRID_GAP_PX = 10

const hashToHue = (value: string) => {
  let hash = 0
  for (let i = 0; i < value.length; i += 1) {
    hash = (hash << 5) - hash + value.charCodeAt(i)
    hash |= 0
  }
  return Math.abs(hash) % 360
}

const paletteFor = (value: string) => {
  const hue = hashToHue(value)
  return {
    bg: `hsl(${hue} 72% 92%)`,
    border: `hsl(${hue} 64% 54%)`,
    text: `hsl(${hue} 45% 28%)`,
    glow: `hsl(${hue} 70% 70% / 0.35)`,
  }
}

const shortLabelFor = (value: string, fallback: string) => {
  const normalized = value.replace(/[^a-z0-9]/gi, '').toUpperCase()
  return normalized.slice(0, 2) || fallback
}

const createGrid = (rows: number, cols: number, fillId: string) =>
  Array.from({ length: rows }, () => Array.from({ length: cols }, () => fillId))

const downloadTextFile = (content: string, filename: string, type: string) => {
  const blob = new Blob([content], { type })
  const url = URL.createObjectURL(blob)
  const link = document.createElement('a')
  link.href = url
  link.download = filename
  link.click()
  URL.revokeObjectURL(url)
}

export function SupercellTool({
  structures,
  onSupercellCreated,
}: SupercellToolProps) {
  const palette = useMemo(() => {
    return structures.map((file, index) => {
      const label = file.label || file.name
      const seed = file.structureId || file.id
      return {
        id: file.structureId ?? null,
        label,
        name: file.name,
        hasLattice: Boolean(file.structure?.lattice),
        palette: paletteFor(seed),
        short: shortLabelFor(label, String.fromCharCode(65 + (index % 26))),
      }
    })
  }, [structures])

  const paletteById = useMemo(() => {
    const map = new Map<string, (typeof palette)[number]>()
    palette.forEach((entry) => {
      if (entry.id) {
        map.set(entry.id, entry)
      }
    })
    return map
  }, [palette])

  const [baseId, setBaseId] = useState<string | null>(null)
  const [brushId, setBrushId] = useState<string | null>(null)
  const [grid, setGrid] = useState<Array<Array<string>>>([])
  const [axisMode, setAxisMode] = useState<AxisMode>('row-b')
  const [checkOverlap, setCheckOverlap] = useState(false)
  const [overlapTolerance, setOverlapTolerance] = useState('0.2')
  const [validateLattice, setValidateLattice] = useState<
    'none' | 'warn' | 'error'
  >('none')
  const [isBuilding, setIsBuilding] = useState(false)
  const [isDownloading, setIsDownloading] = useState(false)
  const [buildError, setBuildError] = useState<string | null>(null)
  const [downloadError, setDownloadError] = useState<string | null>(null)
  const [previewError, setPreviewError] = useState<string | null>(null)
  const [previewStructureId, setPreviewStructureId] = useState<string | null>(
    null,
  )
  const [previewMeta, setPreviewMeta] = useState<SupercellBuildMeta | null>(
    null,
  )

  const baseEntry = baseId ? paletteById.get(baseId) : null
  const baseHasLattice = baseEntry?.hasLattice ?? false
  const gridRows = grid.length
  const gridCols = gridRows > 0 ? grid[0].length : 0

  const gridReady = baseId && gridRows > 0 && gridCols > 0
  const axis: SupercellGridAxis =
    axisMode === 'row-a' ? { row: 'a', col: 'b' } : { row: 'b', col: 'a' }

  const previewUrl = previewStructureId
    ? structureViewUrl(previewStructureId, {
        format: 'cif',
      })
    : null

  const resetGridWithBase = (nextBaseId: string) => {
    const rows = gridRows || DEFAULT_GRID_SIZE
    const cols = gridCols || DEFAULT_GRID_SIZE
    setGrid(createGrid(rows, cols, nextBaseId))
  }

  const handleSelectBase = (nextBaseId: string) => {
    if (!nextBaseId) {
      return
    }
    setBaseId(nextBaseId)
    setBrushId(nextBaseId)
    resetGridWithBase(nextBaseId)
    setBuildError(null)
  }

  const resizeGrid = (nextRows: number, nextCols: number) => {
    if (!baseId) {
      return
    }
    const rows = Math.max(1, nextRows)
    const cols = Math.max(1, nextCols)
    setGrid((prev) =>
      Array.from({ length: rows }, (_row, rowIndex) => {
        const prevRow = prev[rowIndex] ?? []
        return Array.from(
          { length: cols },
          (_col, colIndex) => prevRow[colIndex] ?? baseId,
        )
      }),
    )
  }

  const handlePaintCell = (rowIndex: number, colIndex: number) => {
    if (!brushId) {
      return
    }
    setGrid((prev) =>
      prev.map((row, r) =>
        row.map((cell, c) => {
          if (r !== rowIndex || c !== colIndex) {
            return cell
          }
          return brushId
        }),
      ),
    )
  }

  const handleBuild = async () => {
    setBuildError(null)
    setDownloadError(null)
    setPreviewError(null)
    setPreviewMeta(null)
    setPreviewStructureId(null)
    if (!baseId) {
      setBuildError('Select a base structure.')
      return
    }
    if (!baseHasLattice) {
      setBuildError('Base structure has no lattice data.')
      return
    }
    if (!gridRows || !gridCols) {
      setBuildError('Grid is empty.')
      return
    }
    if (grid.some((row) => row.length !== gridCols)) {
      setBuildError('Grid dimensions are invalid.')
      return
    }
    if (grid.some((row) => row.some((cell) => !cell))) {
      setBuildError('Grid has empty cells.')
      return
    }
    let tolerance: number | undefined
    if (checkOverlap) {
      const parsed = Number(overlapTolerance)
      if (!Number.isFinite(parsed) || parsed < 0) {
        setBuildError('Overlap tolerance is invalid.')
        return
      }
      tolerance = parsed
    }

    setIsBuilding(true)
    try {
      const result = await buildSupercell({
        baseStructureId: baseId,
        grid: {
          rows: gridRows,
          cols: gridCols,
          tiles: grid,
          axis,
        },
        options: {
          checkOverlap,
          overlapTolerance: tolerance,
          validateLattice,
        },
        output: {
          includeStructure: false,
        },
      })
      setPreviewStructureId(result.structureId)
      setPreviewMeta(result.meta)
      setPreviewError(null)
      onSupercellCreated?.({
        structureId: result.structureId,
        meta: result.meta,
      })
    } catch (error) {
      const message =
        error instanceof Error && error.message
          ? error.message
          : 'Failed to build supercell.'
      setBuildError(message)
    } finally {
      setIsBuilding(false)
    }
  }

  const handleDownload = useCallback(async () => {
    setDownloadError(null)
    if (!previewStructureId) {
      setDownloadError('Build a supercell before downloading.')
      return
    }
    setIsDownloading(true)
    try {
      const structure = await getStructure(previewStructureId)
      const xyz = atomsToXyz(structure.atoms)
      const rows = previewMeta?.rows ?? gridRows
      const cols = previewMeta?.cols ?? gridCols
      const timestamp = new Date().toISOString().replace(/[:.]/g, '-')
      const filename = `supercell-${rows}x${cols}-${timestamp}.xyz`
      const atomCount = structure.atoms.length
      const header = `${atomCount}\n${filename}\n`
      downloadTextFile(`${header}${xyz}`, filename, 'chemical/x-xyz')
    } catch (error) {
      const message =
        error instanceof Error && error.message
          ? error.message
          : 'Failed to download XYZ.'
      setDownloadError(message)
    } finally {
      setIsDownloading(false)
    }
  }, [gridCols, gridRows, previewMeta, previewStructureId])

  return (
    <div className="flex h-full gap-4 overflow-hidden">
      <div className="flex w-[26rem] flex-shrink-0 flex-col gap-4 overflow-y-auto pr-1">
        <div className="relative flex aspect-square w-full flex-col overflow-hidden rounded-lg border border-slate-200 bg-white shadow-sm">
          {previewUrl ? (
            <MolstarViewer
              cifUrl={previewUrl}
              onError={setPreviewError}
              onLoad={() => setPreviewError(null)}
            />
          ) : (
            <div className="flex h-full w-full flex-col items-center justify-center text-slate-500">
              <div className="rounded-full bg-indigo-50 p-4">
                <Sparkles className="h-7 w-7 text-indigo-400" />
              </div>
              <p className="mt-3 text-sm font-medium">Supercell Preview</p>
              <p className="text-xs text-slate-400">
                Build to render the result
              </p>
            </div>
          )}
          {previewError ? (
            <div className="absolute inset-0 z-10 flex items-center justify-center bg-white/80 px-4 text-center">
              <div className="rounded-md border border-red-200 bg-red-50 px-3 py-2 text-xs text-red-700 shadow-sm">
                Viewer failed: {previewError}
              </div>
            </div>
          ) : null}
        </div>

        <div className="rounded-lg border border-slate-200 bg-white p-3 text-xs text-slate-600 shadow-sm">
          <div className="flex items-center justify-between gap-2">
            <span className="font-semibold text-slate-700">Build Summary</span>
            <button
              type="button"
              onClick={handleDownload}
              disabled={!previewStructureId || isDownloading}
              className="flex items-center gap-1 rounded-full border border-slate-200 bg-white px-2 py-0.5 text-[10px] font-semibold text-slate-600 shadow-sm transition-colors hover:border-slate-300 hover:text-slate-800 disabled:cursor-not-allowed disabled:border-slate-100 disabled:bg-slate-100 disabled:text-slate-400"
            >
              <Download className="h-3 w-3" />
              {isDownloading ? 'Downloading' : 'Download .xyz'}
            </button>
          </div>
          <div className="mt-3 grid grid-cols-2 gap-2 text-[11px]">
            <div className="rounded-md bg-slate-50 px-2 py-1">
              <div className="text-[10px] uppercase text-slate-400">Size</div>
              <div className="font-semibold text-slate-700">
                {gridRows} x {gridCols}
              </div>
            </div>
            <div className="rounded-md bg-slate-50 px-2 py-1">
              <div className="text-[10px] uppercase text-slate-400">
                Tile Count
              </div>
              <div className="font-semibold text-slate-700">
                {gridRows * gridCols || 0}
              </div>
            </div>
            <div className="rounded-md bg-slate-50 px-2 py-1">
              <div className="text-[10px] uppercase text-slate-400">Axis</div>
              <div className="font-semibold text-slate-700">
                row→{axis.row} / col→{axis.col}
              </div>
            </div>
            <div className="rounded-md bg-slate-50 px-2 py-1">
              <div className="text-[10px] uppercase text-slate-400">
                Overlap
              </div>
              <div className="font-semibold text-slate-700">
                {previewMeta?.overlapCount ?? (checkOverlap ? '...' : '—')}
              </div>
            </div>
          </div>
          {previewMeta ? (
            <div className="mt-3 rounded-md border border-slate-200 bg-white px-2 py-2 text-[10px] text-slate-500">
              base: {previewMeta.baseStructureId.slice(0, 8)} • used:{' '}
              {previewMeta.structureIdsUsed.length}
            </div>
          ) : null}
          {downloadError ? (
            <div className="mt-2 rounded-md border border-red-200 bg-red-50 px-2 py-1 text-[10px] text-red-700">
              {downloadError}
            </div>
          ) : null}
        </div>
      </div>

      <div className="flex min-w-0 flex-1 flex-col overflow-y-auto">
        <div className="flex h-full flex-col overflow-hidden rounded-md border border-blue-200 bg-white shadow-sm">
          <div className="flex shrink-0 items-center gap-2 border-b border-blue-100 bg-blue-50/50 px-3 py-2">
            <Sparkles className="h-4 w-4 text-blue-600" />
            <span className="text-sm font-medium text-blue-900">
              Supercell Builder
            </span>
          </div>
          <div className="flex-1 space-y-5 overflow-y-auto p-4">
            <section className="space-y-2">
              <div className="flex items-center justify-between text-xs font-semibold text-slate-600">
                <span>1. Base structure</span>
                {baseHasLattice ? (
                  <span className="inline-flex items-center gap-1 rounded-full bg-emerald-50 px-2 py-0.5 text-[10px] text-emerald-700">
                    <Check className="h-3 w-3" />
                    lattice ok
                  </span>
                ) : baseId ? (
                  <span className="inline-flex items-center gap-1 rounded-full bg-amber-50 px-2 py-0.5 text-[10px] text-amber-700">
                    <AlertTriangle className="h-3 w-3" />
                    lattice missing
                  </span>
                ) : null}
              </div>
              <div className="relative">
                <select
                  value={baseId ?? ''}
                  onChange={(event) => handleSelectBase(event.target.value)}
                  className="w-full appearance-none rounded-md border border-slate-200 bg-white px-3 py-2 text-sm text-slate-700 shadow-sm focus:border-blue-300 focus:outline-none focus:ring-2 focus:ring-blue-100"
                >
                  <option value="" disabled>
                    Select base structure
                  </option>
                  {palette.map((entry, index) => (
                    <option
                      key={entry.id ?? `label-${entry.label}-${index}`}
                      value={entry.id ?? ''}
                      disabled={!entry.id}
                    >
                      {entry.label}
                      {!entry.id ? ' (no id)' : ''}
                    </option>
                  ))}
                </select>
                <ChevronDown className="pointer-events-none absolute right-3 top-1/2 h-4 w-4 -translate-y-1/2 text-slate-400" />
              </div>
              {!baseId ? (
                <p className="text-[11px] text-slate-400">
                  Base structure must include lattice data.
                </p>
              ) : null}
            </section>

            <section className="space-y-3 border-t border-slate-100 pt-4">
              <div className="flex items-center justify-between text-xs font-semibold text-slate-600">
                <span>2. Grid size &amp; pattern</span>
                <span className="rounded-full bg-slate-100 px-2 py-0.5 text-[10px] text-slate-500">
                  {gridRows} x {gridCols}
                </span>
              </div>
              <GridCanvas
                grid={grid}
                paletteById={paletteById}
                onCellClick={handlePaintCell}
                onResize={resizeGrid}
                disabled={!baseId}
              />
              <p className="text-[11px] text-slate-400">
                Click or drag the dots between tiles to resize. Drag the
                bottom-right handle to expand or shrink.
              </p>
            </section>

            <section className="space-y-3 border-t border-slate-100 pt-4">
              <div className="flex items-center justify-between text-xs font-semibold text-slate-600">
                <span>3. Paint tile</span>
                <span className="text-[10px] text-slate-400">
                  click to replace
                </span>
              </div>
              <div className="grid grid-cols-2 gap-2 sm:grid-cols-3">
                {palette.map((entry, index) => {
                  const disabled = !entry.id
                  const isSelected = entry.id && brushId === entry.id
                  return (
                    <button
                      type="button"
                      key={entry.id ?? `palette-${entry.label}-${index}`}
                      disabled={disabled}
                      onClick={() => entry.id && setBrushId(entry.id)}
                      className={cn(
                        'flex items-center gap-2 rounded-md border px-2 py-2 text-left text-xs font-medium transition',
                        disabled
                          ? 'cursor-not-allowed border-slate-200 bg-slate-50 text-slate-300'
                          : isSelected
                            ? 'border-blue-400 bg-blue-50 text-blue-700 shadow-sm'
                            : 'border-slate-200 bg-white text-slate-600 hover:border-blue-200 hover:text-blue-600',
                      )}
                    >
                      <span
                        className={cn(
                          'flex h-7 w-7 items-center justify-center rounded-md text-[11px] font-semibold',
                        )}
                        style={{
                          backgroundColor: entry.palette.bg,
                          border: `1px solid ${entry.palette.border}`,
                          color: entry.palette.text,
                          boxShadow: `0 0 0 1px ${entry.palette.glow}`,
                        }}
                      >
                        {entry.short}
                      </span>
                      <span className="min-w-0 flex-1 truncate">
                        {entry.label}
                      </span>
                    </button>
                  )
                })}
              </div>
            </section>

            <section className="space-y-3 border-t border-slate-100 pt-4">
              <div className="flex items-center justify-between text-xs font-semibold text-slate-600">
                <span>4. Options</span>
              </div>
              <div className="grid gap-3 sm:grid-cols-2">
                <div className="rounded-md border border-slate-200 bg-slate-50 p-3">
                  <p className="text-[10px] uppercase text-slate-400">
                    Axis mapping
                  </p>
                  <div className="mt-2 grid grid-cols-2 gap-2 text-[11px]">
                    <button
                      type="button"
                      onClick={() => setAxisMode('row-b')}
                      className={cn(
                        'rounded-md border px-2 py-2 text-center font-semibold transition',
                        axisMode === 'row-b'
                          ? 'border-indigo-400 bg-indigo-50 text-indigo-700 shadow-sm'
                          : 'border-slate-200 bg-white text-slate-500 hover:border-indigo-200 hover:text-indigo-600',
                      )}
                    >
                      row→b / col→a
                    </button>
                    <button
                      type="button"
                      onClick={() => setAxisMode('row-a')}
                      className={cn(
                        'rounded-md border px-2 py-2 text-center font-semibold transition',
                        axisMode === 'row-a'
                          ? 'border-indigo-400 bg-indigo-50 text-indigo-700 shadow-sm'
                          : 'border-slate-200 bg-white text-slate-500 hover:border-indigo-200 hover:text-indigo-600',
                      )}
                    >
                      row→a / col→b
                    </button>
                  </div>
                </div>
                <div className="rounded-md border border-slate-200 bg-slate-50 p-3">
                  <p className="text-[10px] uppercase text-slate-400">
                    Validate lattice
                  </p>
                  <div className="mt-2">
                    <select
                      value={validateLattice}
                      onChange={(event) =>
                        setValidateLattice(
                          event.target.value as 'none' | 'warn' | 'error',
                        )
                      }
                      className="w-full rounded-md border border-slate-200 bg-white px-2 py-1.5 text-xs text-slate-600"
                    >
                      <option value="none">none (skip)</option>
                      <option value="warn">warn</option>
                      <option value="error">error</option>
                    </select>
                  </div>
                </div>
                <label className="flex items-center gap-2 rounded-md border border-slate-200 bg-white px-3 py-2 text-xs text-slate-600">
                  <input
                    type="checkbox"
                    checked={checkOverlap}
                    onChange={(event) => setCheckOverlap(event.target.checked)}
                  />
                  overlap check
                </label>
                <div className="rounded-md border border-slate-200 bg-white px-3 py-2 text-xs text-slate-600">
                  <div className="flex items-center justify-between">
                    <span>tolerance</span>
                    <input
                      type="number"
                      min={0}
                      step="0.01"
                      value={overlapTolerance}
                      onChange={(event) =>
                        setOverlapTolerance(event.target.value)
                      }
                      disabled={!checkOverlap}
                      className="w-20 rounded border border-slate-200 px-2 py-1 text-right text-xs text-slate-600 disabled:bg-slate-100"
                    />
                  </div>
                </div>
              </div>
            </section>

            {buildError ? (
              <div className="rounded-md border border-red-200 bg-red-50 px-3 py-2 text-xs text-red-700">
                {buildError}
              </div>
            ) : null}

            <button
              type="button"
              onClick={handleBuild}
              disabled={!gridReady || !baseHasLattice || isBuilding}
              className="flex w-full items-center justify-center gap-2 rounded-md bg-indigo-600 py-2 text-sm font-semibold text-white shadow-sm transition-colors hover:bg-indigo-700 disabled:cursor-not-allowed disabled:bg-slate-300"
            >
              <Sparkles className="h-4 w-4" />
              {isBuilding ? 'Building...' : 'Build Supercell'}
            </button>
          </div>
        </div>
      </div>
    </div>
  )
}

interface GridCanvasProps {
  grid: Array<Array<string>>
  paletteById: Map<
    string,
    { palette: ReturnType<typeof paletteFor>; short: string; label: string }
  >
  onCellClick: (rowIndex: number, colIndex: number) => void
  onResize: (rows: number, cols: number) => void
  disabled?: boolean
}

function GridCanvas({
  grid,
  paletteById,
  onCellClick,
  onResize,
  disabled = false,
}: GridCanvasProps) {
  const rows = grid.length
  const cols = rows > 0 ? grid[0].length : 0
  const renderRows = rows > 0 ? rows + 1 : 0
  const renderCols = cols > 0 ? cols + 1 : 0
  const cell = 'var(--cell-size)'
  const gap = 'var(--gap-size)'
  const rowTracks = Array.from({ length: renderRows * 2 + 1 }, (_, index) =>
    index % 2 === 0 ? gap : cell,
  )
  const colTracks = Array.from({ length: renderCols * 2 + 1 }, (_, index) =>
    index % 2 === 0 ? gap : cell,
  )

  const containerRef = useRef<HTMLDivElement | null>(null)
  const dragRef = useRef<{
    allowBeyond: boolean
    maxRows: number
    maxCols: number
    lastRow: number
    lastCol: number
  } | null>(null)

  const updateFromPointer = useCallback(
    (clientX: number, clientY: number) => {
      if (!containerRef.current) {
        return
      }
      const config = dragRef.current
      if (!config) {
        return
      }
      const container = containerRef.current
      const rect = container.getBoundingClientRect()
      const style = getComputedStyle(container)
      const paddingLeft = Number.parseFloat(style.paddingLeft) || 0
      const paddingTop = Number.parseFloat(style.paddingTop) || 0
      const x = clientX - rect.left - paddingLeft
      const y = clientY - rect.top - paddingTop
      const unit = GRID_CELL_PX + GRID_GAP_PX
      const base = GRID_GAP_PX / 2
      let nextRow = Math.round((y - base) / unit)
      let nextCol = Math.round((x - base) / unit)
      nextRow = Math.max(1, nextRow)
      nextCol = Math.max(1, nextCol)
      if (!config.allowBeyond) {
        nextRow = Math.min(nextRow, config.maxRows)
        nextCol = Math.min(nextCol, config.maxCols)
      }
      if (config.lastRow === nextRow && config.lastCol === nextCol) {
        return
      }
      dragRef.current = {
        ...config,
        lastRow: nextRow,
        lastCol: nextCol,
      }
      onResize(nextRow, nextCol)
    },
    [onResize],
  )

  useEffect(() => {
    const handleMove = (event: PointerEvent) => {
      if (!dragRef.current) {
        return
      }
      updateFromPointer(event.clientX, event.clientY)
    }
    const handleUp = () => {
      dragRef.current = null
    }
    window.addEventListener('pointermove', handleMove)
    window.addEventListener('pointerup', handleUp)
    window.addEventListener('pointercancel', handleUp)
    return () => {
      window.removeEventListener('pointermove', handleMove)
      window.removeEventListener('pointerup', handleUp)
      window.removeEventListener('pointercancel', handleUp)
    }
  }, [updateFromPointer])

  const handleDotPointerDown = (
    event: ReactPointerEvent<HTMLButtonElement>,
    row: number,
    col: number,
    allowBeyond: boolean,
  ) => {
    if (disabled) {
      return
    }
    event.preventDefault()
    dragRef.current = {
      allowBeyond,
      maxRows: renderRows,
      maxCols: renderCols,
      lastRow: row,
      lastCol: col,
    }
    onResize(row, col)
    updateFromPointer(event.clientX, event.clientY)
  }

  if (!rows || !cols) {
    return (
      <div className="rounded-lg border border-slate-200 bg-slate-50 p-6 text-center text-[11px] text-slate-400">
        Select a base structure to edit the grid.
      </div>
    )
  }

  return (
    <div
      ref={containerRef}
      className={cn(
        'relative grid rounded-lg border border-slate-200 bg-slate-50 p-2',
        disabled && 'opacity-60',
      )}
      style={
        {
          gridTemplateRows: rowTracks.join(' '),
          gridTemplateColumns: colTracks.join(' '),
          '--cell-size': `${GRID_CELL_PX}px`,
          '--gap-size': `${GRID_GAP_PX}px`,
        } as CSSProperties
      }
    >
      {Array.from({ length: renderRows }).map((_, rowIndex) =>
        Array.from({ length: renderCols }).map((__, colIndex) => {
          const isGhost = rowIndex >= rows || colIndex >= cols
          const cellId = !isGhost ? grid[rowIndex]?.[colIndex] : null
          const entry = cellId ? paletteById.get(cellId) : null
          const palette = cellId ? (entry?.palette ?? paletteFor(cellId)) : null
          return (
            <button
              key={`cell-${rowIndex}-${colIndex}`}
              type="button"
              disabled={disabled || isGhost}
              onClick={() => {
                if (!isGhost) {
                  onCellClick(rowIndex, colIndex)
                }
              }}
              className={cn(
                'relative flex items-center justify-center rounded-md border text-[11px] font-semibold uppercase transition',
                isGhost
                  ? 'cursor-default border-dashed border-slate-200 bg-white/60 text-slate-300'
                  : disabled
                    ? 'cursor-not-allowed border-slate-200 bg-white/70 text-slate-300'
                    : 'border-slate-200 hover:shadow-sm',
              )}
              style={{
                gridRow: rowIndex * 2 + 2,
                gridColumn: colIndex * 2 + 2,
                backgroundColor: isGhost
                  ? 'rgba(248,250,252,0.8)'
                  : palette?.bg,
                borderColor: isGhost ? '#e2e8f0' : palette?.border,
                color: palette?.text,
                boxShadow: isGhost ? 'none' : `0 0 0 1px ${palette?.glow}`,
              }}
              title={
                isGhost ? 'Expansion preview' : (entry?.label ?? cellId ?? '')
              }
            >
              {isGhost ? '' : (entry?.short ?? '•')}
            </button>
          )
        }),
      )}

      {Array.from({ length: renderRows }).map((_, boundaryRow) =>
        Array.from({ length: renderCols }).map((__, boundaryCol) => {
          const row = boundaryRow + 1
          const col = boundaryCol + 1
          const isHandle = row === renderRows && col === renderCols
          return (
            <DotControl
              key={`dot-${row}-${col}`}
              row={row}
              col={col}
              gridRow={row * 2 + 1}
              gridColumn={col * 2 + 1}
              onPointerDown={handleDotPointerDown}
              disabled={disabled}
              isHandle={isHandle}
            />
          )
        }),
      )}
    </div>
  )
}

interface DotControlProps {
  row: number
  col: number
  gridRow: number | string
  gridColumn: number | string
  disabled: boolean
  isHandle?: boolean
  onPointerDown: (
    event: ReactPointerEvent<HTMLButtonElement>,
    row: number,
    col: number,
    allowBeyond: boolean,
  ) => void
}

function DotControl({
  row,
  col,
  gridRow,
  gridColumn,
  disabled,
  isHandle = false,
  onPointerDown,
}: DotControlProps) {
  return (
    <button
      type="button"
      disabled={disabled}
      onPointerDown={(event) =>
        onPointerDown(event, row, col, Boolean(isHandle))
      }
      className={cn(
        'group z-10 flex items-center justify-center rounded-full border transition',
        disabled
          ? 'cursor-not-allowed border-slate-200 bg-white/70'
          : isHandle
            ? 'border-indigo-300 bg-white text-indigo-500 shadow-sm hover:border-indigo-400 hover:text-indigo-600'
            : 'border-slate-200 bg-white text-slate-400 hover:border-slate-300 hover:text-slate-500',
      )}
      style={{
        gridRow,
        gridColumn,
        width: '0.9rem',
        height: '0.9rem',
        justifySelf: 'center',
        alignSelf: 'center',
      }}
      title={`${row} x ${col}`}
    >
      {isHandle ? (
        <MoveDiagonal2 className="h-3 w-3" />
      ) : (
        <span className="h-1.5 w-1.5 rounded-full bg-current" />
      )}
    </button>
  )
}
