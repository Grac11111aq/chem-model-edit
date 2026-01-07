import { useEffect, useMemo, useRef, useState } from 'react'

import { AtomEditorPanel } from './components/AtomEditorPanel'
import { ComparePanel } from './components/ComparePanel'
import { ExportSettings } from './components/ExportSettings'
import { SelectionTools } from './components/SelectionTools'
import { StructureSidebar } from './components/StructureSidebar'
import { ViewControls } from './components/ViewControls'
import { ViewerPanel } from './components/ViewerPanel'
import {
  LATTICE_TOLERANCE,
  invert3,
  latticeToMatrix,
  maxAbsDifference,
  minimumImageDelta,
} from './lattice'
import type { MolstarStructure } from './components/types'

import type { Atom } from '@/lib/types'
import type { AtomDraft, PbcState, StructureState } from './types'
import { atomsToXyz, parseXyzBlock } from '@/lib/xyz'
import { atomsToPdb } from '@/lib/pdb'
import { exportQeInput, parseQeInput } from '@/lib/api'
import {
  alignSelectedCentroid,
  alignSelectedToOrigin,
  shiftAtoms,
} from '@/components/compare/align'
import { downloadShareHtml } from '@/components/share/html-export'

const palette = [
  'from-sky-300 to-cyan-300',
  'from-rose-300 to-orange-300',
  'from-lime-300 to-emerald-300',
  'from-amber-300 to-yellow-200',
  'from-fuchsia-300 to-purple-300',
]

const sampleQe = [
  '&CONTROL',
  "  calculation='scf'",
  '/',
  '&SYSTEM',
  '  ibrav=0, nat=3, ntyp=2',
  '/',
  '&ELECTRONS',
  '/',
  'ATOMIC_SPECIES',
  ' O 15.999 O.pbe-rrkjus.UPF',
  ' H 1.0079 H.pbe-rrkjus.UPF',
  'CELL_PARAMETERS angstrom',
  '  6.0 0.0 0.0',
  '  0.0 6.0 0.0',
  '  0.0 0.0 6.0',
  'ATOMIC_POSITIONS angstrom',
  ' O 0.0000 0.0000 0.0000',
  ' H 0.7570 0.5860 0.0000',
  ' H -0.7570 0.5860 0.0000',
].join('\n')

const samplePdb = [
  'HEADER    SAMPLE',
  'ATOM      1  O   HOH A   1       0.000   0.000   0.000  1.00  0.00           O',
  'ATOM      2  H1  HOH A   1       0.757   0.586   0.000  1.00  0.00           H',
  'ATOM      3  H2  HOH A   1      -0.757   0.586   0.000  1.00  0.00           H',
  'END',
].join('\n')

const makeDrafts = (nextAtoms: Array<Atom>) =>
  nextAtoms.map((atom) => ({
    symbol: atom.symbol,
    x: atom.x.toFixed(4),
    y: atom.y.toFixed(4),
    z: atom.z.toFixed(4),
  }))

export default function EditorPage() {
  const [structures, setStructures] = useState<Array<StructureState>>([
    {
      id: 'A',
      name: 'Catalyst-A',
      color: palette[0],
      atoms: [],
      drafts: [],
      isVisible: true,
      opacity: 1,
      lattice: null,
    },
    {
      id: 'B',
      name: 'Ligand-B',
      color: palette[1],
      atoms: [],
      drafts: [],
      isVisible: true,
      opacity: 0.5,
      lattice: null,
    },
  ])
  const [activeId, setActiveId] = useState(structures[0]?.id ?? 'A')
  const [overlayEnabled, setOverlayEnabled] = useState(true)
  const [qeInput, setQeInput] = useState(sampleQe)
  const [selectedIndices, setSelectedIndices] = useState<Array<number>>([])
  const [shiftDraft, setShiftDraft] = useState({
    x: '0.000',
    y: '0.000',
    z: '0.000',
  })
  const [isParsing, setIsParsing] = useState(false)
  const [error, setError] = useState<string | null>(null)
  const [exported, setExported] = useState('')
  const fileInputRef = useRef<HTMLInputElement | null>(null)
  const [compareTargetId, setCompareTargetId] = useState<string | null>(null)

  const activeStructure =
    structures.find((structure) => structure.id === activeId) ?? structures[0]
  const atoms = activeStructure.atoms
  const atomDrafts = activeStructure.drafts

  const overlayTargets = useMemo(() => {
    if (overlayEnabled) {
      return structures
    }
    return structures.filter((structure) => structure.id === activeId)
  }, [activeId, overlayEnabled, structures])

  const viewerStructures = useMemo<Array<MolstarStructure>>(() => {
    return overlayTargets.flatMap((structure) => {
      if (!structure.isVisible) {
        return []
      }
      const pdbText =
        structure.atoms.length > 0
          ? atomsToPdb(structure.atoms)
          : structure.id === activeId
            ? samplePdb
            : null
      if (!pdbText) {
        return []
      }
      return [
        {
          id: structure.id,
          pdbText,
          opacity: structure.opacity,
          visible: structure.isVisible,
        },
      ]
    })
  }, [activeId, overlayTargets])

  const visibleOverlayCount = useMemo(
    () => overlayTargets.filter((structure) => structure.isVisible).length,
    [overlayTargets],
  )

  const compareCandidates = useMemo(
    () => structures.filter((structure) => structure.id !== activeId),
    [activeId, structures],
  )

  const compareTarget = useMemo(
    () =>
      structures.find((structure) => structure.id === compareTargetId) ?? null,
    [compareTargetId, structures],
  )

  const countMismatch = useMemo(() => {
    if (!compareTarget) {
      return false
    }
    if (atoms.length === 0 || compareTarget.atoms.length === 0) {
      return false
    }
    return atoms.length !== compareTarget.atoms.length
  }, [atoms.length, compareTarget])

  const pbcState = useMemo<PbcState>(() => {
    if (!compareTarget) {
      return { enabled: false, reason: '比較対象なし' }
    }
    const latticeA = activeStructure.lattice ?? null
    const latticeB = compareTarget.lattice ?? null
    if (!latticeA || !latticeB) {
      return { enabled: false, reason: '格子情報が不足' }
    }
    const matrixA = latticeToMatrix(latticeA)
    const matrixB = latticeToMatrix(latticeB)
    const mismatch = maxAbsDifference(matrixA, matrixB) > LATTICE_TOLERANCE
    if (mismatch) {
      return { enabled: false, reason: '格子が一致しない', mismatch: true }
    }
    const inverse = invert3(matrixA)
    if (!inverse) {
      return { enabled: false, reason: '格子が特異', mismatch: false }
    }
    return { enabled: true, reason: 'PBC有効', matrix: matrixA, inverse }
  }, [activeStructure.lattice, compareTarget])

  const selectedAtoms = useMemo(
    () =>
      selectedIndices
        .filter((index) => index >= 0 && index < atoms.length)
        .map((index) => atoms[index]),
    [atoms, selectedIndices],
  )

  const distanceRows = useMemo(() => {
    if (
      !compareTarget ||
      atoms.length === 0 ||
      compareTarget.atoms.length === 0
    ) {
      return []
    }
    const limit = Math.min(atoms.length, compareTarget.atoms.length)
    const indices =
      selectedIndices.length > 0
        ? selectedIndices.filter((index) => index >= 0 && index < limit)
        : [...Array(Math.min(limit, 200)).keys()]
    return indices.map((index) => {
      const source = atoms[index]
      const target = compareTarget.atoms[index]
      const delta = {
        x: source.x - target.x,
        y: source.y - target.y,
        z: source.z - target.z,
      }
      const adjusted = pbcState.enabled
        ? minimumImageDelta(delta, pbcState.matrix, pbcState.inverse)
        : delta
      return {
        index,
        source,
        target,
        dist: Math.sqrt(adjusted.x ** 2 + adjusted.y ** 2 + adjusted.z ** 2),
      }
    })
  }, [atoms, compareTarget, selectedIndices, pbcState])

  const selectedDistance = useMemo(() => {
    if (selectedAtoms.length < 2) {
      return null
    }
    const [a, b] = selectedAtoms
    const dx = a.x - b.x
    const dy = a.y - b.y
    const dz = a.z - b.z
    return Math.sqrt(dx * dx + dy * dy + dz * dz)
  }, [selectedAtoms])

  const updateActive = (
    updater: (structure: StructureState) => StructureState,
  ) => {
    setStructures((prev) =>
      prev.map((structure) =>
        structure.id === activeId ? updater(structure) : structure,
      ),
    )
  }

  const updateStructure = (
    id: string,
    updater: (structure: StructureState) => StructureState,
  ) => {
    setStructures((prev) =>
      prev.map((structure) =>
        structure.id === id ? updater(structure) : structure,
      ),
    )
  }

  const syncDrafts = (nextAtoms: Array<Atom>) => {
    updateActive((structure) => ({
      ...structure,
      drafts: makeDrafts(nextAtoms),
    }))
  }

  const parseContent = async (content: string) => {
    setIsParsing(true)
    setError(null)
    setExported('')
    try {
      const structure = await parseQeInput(content)
      updateActive((current) => ({
        ...current,
        atoms: structure.atoms,
        drafts: makeDrafts(structure.atoms),
        lattice: structure.lattice ?? null,
      }))
      setSelectedIndices([])
    } catch (err) {
      setError(err instanceof Error ? err.message : 'パースに失敗しました。')
    } finally {
      setIsParsing(false)
    }
  }

  const handleParse = async () => {
    await parseContent(qeInput)
  }

  const handleExport = async () => {
    if (atoms.length === 0) {
      setError('エクスポートする原子がありません。')
      return
    }
    setError(null)
    try {
      const content = await exportQeInput({ atoms })
      setExported(content)
      await navigator.clipboard.writeText(content)
    } catch (err) {
      setError(
        err instanceof Error ? err.message : 'エクスポートに失敗しました。',
      )
    }
  }

  const handleShareHtml = () => {
    const shareTargets = overlayTargets.filter(
      (structure) => structure.isVisible && structure.atoms.length > 0,
    )
    if (shareTargets.length === 0) {
      setError('共有する原子がありません。')
      return
    }
    setError(null)
    downloadShareHtml(structures, { activeId, overlayEnabled })
  }

  const handleAtomChange = (
    index: number,
    field: keyof AtomDraft,
    value: string,
  ) => {
    updateActive((structure) => {
      const nextDrafts = [...structure.drafts]
      const currentDraft = nextDrafts[index] ?? {
        symbol: '',
        x: '',
        y: '',
        z: '',
      }
      nextDrafts[index] = { ...currentDraft, [field]: value }

      const nextAtoms = [...structure.atoms]
      if (index < 0 || index >= nextAtoms.length) {
        return { ...structure, drafts: nextDrafts }
      }

      const currentAtom = nextAtoms[index]
      if (field === 'symbol') {
        const trimmed = value.trim()
        nextAtoms[index] = {
          ...currentAtom,
          symbol: trimmed || currentAtom.symbol,
        }
        return { ...structure, drafts: nextDrafts, atoms: nextAtoms }
      }

      const parsed = Number(value)
      if (Number.isNaN(parsed)) {
        return { ...structure, drafts: nextDrafts }
      }

      nextAtoms[index] = { ...currentAtom, [field]: parsed }
      return { ...structure, drafts: nextDrafts, atoms: nextAtoms }
    })
  }

  const handleAddAtom = () => {
    const nextAtom: Atom = { symbol: 'X', x: 0, y: 0, z: 0 }
    updateActive((structure) => ({
      ...structure,
      atoms: [...structure.atoms, nextAtom],
      drafts: [
        ...structure.drafts,
        { symbol: 'X', x: '0.0000', y: '0.0000', z: '0.0000' },
      ],
    }))
  }

  const handleToggleVisibility = (id: string, next: boolean) => {
    updateStructure(id, (structure) => ({ ...structure, isVisible: next }))
  }

  const handleOpacityChange = (id: string, value: number) => {
    const nextOpacity = Math.min(1, Math.max(0, value))
    updateStructure(id, (structure) => ({ ...structure, opacity: nextOpacity }))
  }

  const toggleSelect = (index: number) => {
    setSelectedIndices((prev) =>
      prev.includes(index)
        ? prev.filter((item) => item !== index)
        : [...prev, index],
    )
  }

  const handleCopySelected = async () => {
    if (selectedAtoms.length === 0) {
      setError('コピーする原子を選択してください。')
      return
    }
    setError(null)
    try {
      await navigator.clipboard.writeText(atomsToXyz(selectedAtoms))
    } catch (err) {
      setError(err instanceof Error ? err.message : 'コピーに失敗しました。')
    }
  }

  const handlePasteAppend = async () => {
    try {
      const text = await navigator.clipboard.readText()
      const parsed = parseXyzBlock(text)
      if (parsed.length === 0) {
        setError('クリップボードに有効なXYZがありません。')
        return
      }
      setError(null)
      const start = atoms.length
      updateActive((structure) => ({
        ...structure,
        atoms: [...structure.atoms, ...parsed],
        drafts: [...structure.drafts, ...makeDrafts(parsed)],
      }))
      setSelectedIndices(parsed.map((_, i) => start + i))
    } catch (err) {
      setError(err instanceof Error ? err.message : '貼り付けに失敗しました。')
    }
  }

  const handleShiftSelected = () => {
    if (selectedIndices.length === 0) {
      setError('シフトする原子を選択してください。')
      return
    }
    const dx = Number(shiftDraft.x)
    const dy = Number(shiftDraft.y)
    const dz = Number(shiftDraft.z)
    if ([dx, dy, dz].some((val) => Number.isNaN(val))) {
      setError('シフト量の数値が無効です。')
      return
    }
    setError(null)
    const nextAtoms = shiftAtoms(atoms, selectedIndices, {
      x: dx,
      y: dy,
      z: dz,
    })
    syncDrafts(nextAtoms)
    updateActive((structure) => ({ ...structure, atoms: nextAtoms }))
  }

  const handleAlignOrigin = () => {
    if (selectedIndices.length === 0) {
      setError('整列する原子を選択してください。')
      return
    }
    setError(null)
    const nextAtoms = alignSelectedToOrigin(atoms, selectedIndices)
    syncDrafts(nextAtoms)
    updateActive((structure) => ({ ...structure, atoms: nextAtoms }))
  }

  const handleAlignCentroid = () => {
    if (selectedIndices.length === 0) {
      setError('整列する原子を選択してください。')
      return
    }
    setError(null)
    const nextAtoms = alignSelectedCentroid(atoms, selectedIndices)
    syncDrafts(nextAtoms)
    updateActive((structure) => ({ ...structure, atoms: nextAtoms }))
  }

  const handleSurfaceTransfer = () => {
    if (!compareTarget) {
      setError('転写先の構造がありません。')
      return
    }
    if (atoms.length === 0 || compareTarget.atoms.length === 0) {
      setError('転写元・転写先の両方に構造を読み込んでください。')
      return
    }
    if (selectedIndices.length === 0) {
      setError('転写する原子を選択してください。')
      return
    }
    const invalid = selectedIndices.filter(
      (index) => index >= compareTarget.atoms.length,
    )
    if (invalid.length > 0) {
      setError('転写先の原子数が不足しています。')
      return
    }
    setError(null)
    updateStructure(compareTarget.id, (structure) => {
      const nextAtoms = structure.atoms.map((atom, index) => {
        if (!selectedIndices.includes(index)) {
          return atom
        }
        const sourceAtom = atoms[index]
        return {
          ...atom,
          x: sourceAtom.x,
          y: sourceAtom.y,
          z: sourceAtom.z,
        }
      })
      return { ...structure, atoms: nextAtoms, drafts: makeDrafts(nextAtoms) }
    })
  }

  const handleSelectStructure = (id: string) => {
    setActiveId(id)
    setSelectedIndices([])
    setError(null)
    setExported('')
  }

  const handleAddStructure = () => {
    setStructures((prev) => {
      const index = prev.length + 1
      const id = index <= 26 ? String.fromCharCode(64 + index) : `S${index}`
      const newStructure: StructureState = {
        id,
        name: `Structure-${id}`,
        color: palette[(index - 1) % palette.length],
        atoms: [],
        drafts: [],
        isVisible: true,
        opacity: 1,
        lattice: null,
      }
      setActiveId(id)
      setSelectedIndices([])
      return [...prev, newStructure]
    })
  }

  const handleImportFile = async (
    event: React.ChangeEvent<HTMLInputElement>,
  ) => {
    const file = event.target.files?.[0]
    if (!file) {
      return
    }
    const text = await file.text()
    setQeInput(text)
    await parseContent(text)
    event.target.value = ''
  }

  const handleImportClipboard = async () => {
    try {
      const text = await navigator.clipboard.readText()
      if (!text.trim()) {
        setError('クリップボードが空です。')
        return
      }
      setQeInput(text)
      await parseContent(text)
    } catch (err) {
      setError(
        err instanceof Error
          ? err.message
          : 'クリップボード読取に失敗しました。',
      )
    }
  }

  useEffect(() => {
    const handler = () => fileInputRef.current?.click()
    window.addEventListener('chem-model-import', handler)
    return () => window.removeEventListener('chem-model-import', handler)
  }, [])

  useEffect(() => {
    if (compareCandidates.length === 0) {
      setCompareTargetId(null)
      return
    }
    setCompareTargetId((prev) => {
      if (
        prev &&
        compareCandidates.some((candidate) => candidate.id === prev)
      ) {
        return prev
      }
      return compareCandidates[0].id
    })
  }, [compareCandidates])

  return (
    <div className="min-h-screen bg-slate-950 text-white">
      <div className="relative overflow-hidden">
        <div className="pointer-events-none absolute inset-0 bg-[radial-gradient(circle_at_top,_rgba(248,186,94,0.16),_transparent_55%)]" />
        <div className="pointer-events-none absolute inset-0 bg-[radial-gradient(circle_at_20%_20%,_rgba(56,189,248,0.08),_transparent_45%)]" />
        <div className="pointer-events-none absolute inset-0 bg-[radial-gradient(circle_at_80%_20%,_rgba(251,113,133,0.12),_transparent_45%)]" />
        <main className="relative mx-auto flex min-h-[calc(100vh-64px)] max-w-[1400px] flex-col gap-6 px-6 py-6 lg:flex-row">
          <section className="flex w-full flex-col gap-4 lg:max-w-[260px]">
            <StructureSidebar
              structures={structures}
              activeId={activeId}
              overlayEnabled={overlayEnabled}
              fileInputRef={fileInputRef}
              onToggleOverlay={() => setOverlayEnabled((prev) => !prev)}
              onAddStructure={handleAddStructure}
              onSelectStructure={handleSelectStructure}
              onToggleVisibility={handleToggleVisibility}
              onOpacityChange={handleOpacityChange}
              onImportFile={handleImportFile}
              onImportClipboard={handleImportClipboard}
            />
            <ComparePanel
              compareCandidates={compareCandidates}
              compareTarget={compareTarget}
              countMismatch={countMismatch}
              pbcState={pbcState}
              distanceRows={distanceRows}
              selectedIndices={selectedIndices}
              atomsCount={atoms.length}
              selectedDistance={selectedDistance}
              activeId={activeId}
              onChangeTarget={setCompareTargetId}
              onAlignSelected={handleAlignCentroid}
              onSurfaceTransfer={handleSurfaceTransfer}
            />
          </section>

          <section className="flex w-full flex-1 flex-col gap-4">
            <AtomEditorPanel
              qeInput={qeInput}
              isParsing={isParsing}
              error={error}
              exported={exported}
              atoms={atoms}
              atomDrafts={atomDrafts}
              selectedIndices={selectedIndices}
              onChangeQeInput={setQeInput}
              onParse={handleParse}
              onResetSample={() => setQeInput(sampleQe)}
              onExport={handleExport}
              onShare={handleShareHtml}
              onAddAtom={handleAddAtom}
              onToggleSelect={toggleSelect}
              onAtomChange={handleAtomChange}
            />

            <div className="grid gap-4 lg:grid-cols-2">
              <SelectionTools
                shiftDraft={shiftDraft}
                onShiftDraftChange={(axis, value) =>
                  setShiftDraft((prev) => ({ ...prev, [axis]: value }))
                }
                onApplyShift={handleShiftSelected}
                onAlignOrigin={handleAlignOrigin}
                onAlignCentroid={handleAlignCentroid}
                onCopySelected={handleCopySelected}
                onPasteAppend={handlePasteAppend}
                onClearSelection={() => setSelectedIndices([])}
              />
              <ViewControls
                overlayEnabled={overlayEnabled}
                visibleOverlayCount={visibleOverlayCount}
                activeOpacity={activeStructure.opacity}
              />
            </div>

            <ExportSettings />
          </section>

          <section className="flex w-full flex-col gap-4 lg:max-w-[360px]">
            <ViewerPanel structures={viewerStructures} />
          </section>
        </main>
      </div>
    </div>
  )
}
