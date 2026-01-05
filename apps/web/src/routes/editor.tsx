import { createFileRoute } from '@tanstack/react-router'
import { Fragment, useEffect, useMemo, useRef, useState } from 'react'
import {
  ArrowDownToLine,
  ChevronDown,
  Droplets,
  Layers,
  PencilRuler,
  Plus,
  RefreshCw,
} from 'lucide-react'
import MolstarViewer from '../components/molstar/MolstarViewer'
import { downloadShareHtml } from '../components/share/html-export'
import { exportQeInput, parseQeInput } from '../lib/api'
import { atomsToPdb } from '../lib/pdb'
import type { Atom, Lattice, Vector3 } from '../lib/types'
import { atomsToXyz, parseXyzBlock } from '../lib/xyz'
import {
  alignSelectedCentroid,
  alignSelectedToOrigin,
  shiftAtoms,
} from '../components/compare/align'

type AtomDraft = {
  symbol: string
  x: string
  y: string
  z: string
}

type StructureState = {
  id: string
  name: string
  color: string
  atoms: Atom[]
  drafts: AtomDraft[]
  isVisible: boolean
  opacity: number
  lattice?: Lattice | null
}

type Matrix3 = [
  [number, number, number],
  [number, number, number],
  [number, number, number],
]

const LATTICE_TOLERANCE = 1.0e-3

const latticeToMatrix = (lattice: Lattice): Matrix3 => [
  [lattice.a.x, lattice.a.y, lattice.a.z],
  [lattice.b.x, lattice.b.y, lattice.b.z],
  [lattice.c.x, lattice.c.y, lattice.c.z],
]

const maxAbsDifference = (a: Matrix3, b: Matrix3) => {
  let max = 0
  for (let i = 0; i < 3; i += 1) {
    for (let j = 0; j < 3; j += 1) {
      const diff = Math.abs(a[i][j] - b[i][j])
      if (diff > max) {
        max = diff
      }
    }
  }
  return max
}

const invert3 = (m: Matrix3): Matrix3 | null => {
  const [
    [a00, a01, a02],
    [a10, a11, a12],
    [a20, a21, a22],
  ] = m
  const b01 = a22 * a11 - a12 * a21
  const b11 = -a22 * a10 + a12 * a20
  const b21 = a21 * a10 - a11 * a20
  const det = a00 * b01 + a01 * b11 + a02 * b21
  if (Math.abs(det) < 1.0e-12) {
    return null
  }
  const invDet = 1 / det
  return [
    [
      b01 * invDet,
      (-a22 * a01 + a02 * a21) * invDet,
      (a12 * a01 - a02 * a11) * invDet,
    ],
    [
      b11 * invDet,
      (a22 * a00 - a02 * a20) * invDet,
      (-a12 * a00 + a02 * a10) * invDet,
    ],
    [
      b21 * invDet,
      (-a21 * a00 + a01 * a20) * invDet,
      (a11 * a00 - a01 * a10) * invDet,
    ],
  ]
}

const multiplyMatVec = (m: Matrix3, v: Vector3): Vector3 => ({
  x: m[0][0] * v.x + m[0][1] * v.y + m[0][2] * v.z,
  y: m[1][0] * v.x + m[1][1] * v.y + m[1][2] * v.z,
  z: m[2][0] * v.x + m[2][1] * v.y + m[2][2] * v.z,
})

const wrapFractional = (v: Vector3): Vector3 => ({
  x: v.x - Math.round(v.x),
  y: v.y - Math.round(v.y),
  z: v.z - Math.round(v.z),
})

const minimumImageDelta = (
  delta: Vector3,
  lattice: Matrix3,
  inverse: Matrix3,
): Vector3 => {
  const fractional = multiplyMatVec(inverse, delta)
  const wrapped = wrapFractional(fractional)
  return multiplyMatVec(lattice, wrapped)
}

export const Route = createFileRoute('/editor')({
  component: EditorPage,
})

function EditorPage() {
  const palette = [
    'from-sky-300 to-cyan-300',
    'from-rose-300 to-orange-300',
    'from-lime-300 to-emerald-300',
    'from-amber-300 to-yellow-200',
    'from-fuchsia-300 to-purple-300',
  ]
  const [structures, setStructures] = useState<StructureState[]>([
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
  const [qeInput, setQeInput] = useState(sampleQe)
  const [selectedIndices, setSelectedIndices] = useState<number[]>([])
  const [shiftDraft, setShiftDraft] = useState({ x: '0.000', y: '0.000', z: '0.000' })
  const [isParsing, setIsParsing] = useState(false)
  const [error, setError] = useState<string | null>(null)
  const [exported, setExported] = useState('')
  const fileInputRef = useRef<HTMLInputElement | null>(null)
  const [compareTargetId, setCompareTargetId] = useState<string | null>(null)

  const activeStructure =
    structures.find((structure) => structure.id === activeId) ?? structures[0]
  const atoms = activeStructure?.atoms ?? []
  const atomDrafts = activeStructure?.drafts ?? []

  const samplePdb = [
    'HEADER    SAMPLE',
    'ATOM      1  O   HOH A   1       0.000   0.000   0.000  1.00  0.00           O',
    'ATOM      2  H1  HOH A   1       0.757   0.586   0.000  1.00  0.00           H',
    'ATOM      3  H2  HOH A   1      -0.757   0.586   0.000  1.00  0.00           H',
    'END',
  ].join('\n')
  const overlayTargets = useMemo(() => {
    if (overlayEnabled) {
      return structures
    }
    return structures.filter((structure) => structure.id === activeId)
  }, [activeId, overlayEnabled, structures])
  const viewerStructures = useMemo(() => {
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
  }, [activeId, overlayTargets, samplePdb])
  const visibleOverlayCount = useMemo(
    () => overlayTargets.filter((structure) => structure.isVisible).length,
    [overlayTargets],
  )
  const compareCandidates = useMemo(
    () => structures.filter((structure) => structure.id !== activeId),
    [activeId, structures],
  )
  const compareTarget = useMemo(
    () => structures.find((structure) => structure.id === compareTargetId) ?? null,
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
  const pbcState = useMemo(() => {
    if (!compareTarget) {
      return { enabled: false, reason: '比較対象なし' }
    }
    const latticeA = activeStructure?.lattice ?? null
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
  }, [activeStructure?.lattice, compareTarget, activeStructure])
  const selectedAtoms = useMemo(
    () => selectedIndices.map((index) => atoms[index]).filter(Boolean),
    [atoms, selectedIndices],
  )
  const distanceRows = useMemo(() => {
    if (!compareTarget || atoms.length === 0 || compareTarget.atoms.length === 0) {
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
      const adjusted =
        pbcState.enabled && pbcState.matrix && pbcState.inverse
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

  const updateActive = (updater: (structure: StructureState) => StructureState) => {
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
      prev.map((structure) => (structure.id === id ? updater(structure) : structure)),
    )
  }

  const makeDrafts = (nextAtoms: Atom[]) =>
    nextAtoms.map((atom) => ({
      symbol: atom.symbol,
      x: atom.x.toFixed(4),
      y: atom.y.toFixed(4),
      z: atom.z.toFixed(4),
    }))

  const syncDrafts = (nextAtoms: Atom[]) => {
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
      setError(err instanceof Error ? err.message : 'エクスポートに失敗しました。')
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
      const current = nextDrafts[index] ?? { symbol: '', x: '', y: '', z: '' }
      nextDrafts[index] = { ...current, [field]: value }
      return { ...structure, drafts: nextDrafts }
    })

    if (field === 'symbol') {
      updateActive((structure) => {
        const nextAtoms = [...structure.atoms]
        const current = nextAtoms[index]
        if (!current) {
          return structure
        }
        nextAtoms[index] = {
          ...current,
          symbol: value.trim() || current.symbol,
        }
        return { ...structure, atoms: nextAtoms }
      })
      return
    }

    const parsed = Number(value)
    if (Number.isNaN(parsed)) {
      return
    }
    updateActive((structure) => {
      const nextAtoms = [...structure.atoms]
      const current = nextAtoms[index]
      if (!current) {
        return structure
      }
      nextAtoms[index] = { ...current, [field]: parsed }
      return { ...structure, atoms: nextAtoms }
    })
  }

  const handleAddAtom = () => {
    const nextAtom: Atom = { symbol: 'X', x: 0, y: 0, z: 0 }
    updateActive((structure) => ({
      ...structure,
      atoms: [...structure.atoms, nextAtom],
      drafts: [...structure.drafts, { symbol: 'X', x: '0.0000', y: '0.0000', z: '0.0000' }],
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
      prev.includes(index) ? prev.filter((item) => item !== index) : [...prev, index],
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
    const nextAtoms = shiftAtoms(atoms, selectedIndices, { x: dx, y: dy, z: dz })
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
    const invalid = selectedIndices.filter((index) => index >= compareTarget.atoms.length)
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
      const id =
        index <= 26 ? String.fromCharCode(64 + index) : `S${index}`
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
      setError(err instanceof Error ? err.message : 'クリップボード読取に失敗しました。')
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
      if (prev && compareCandidates.some((candidate) => candidate.id === prev)) {
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
            <div className="rounded-2xl border border-white/10 bg-white/5 p-4 shadow-lg shadow-black/30">
              <div className="flex items-center justify-between">
                <p className="text-xs uppercase tracking-[0.3em] text-white/50">
                  Structures
                </p>
                <div className="flex items-center gap-2">
                  <button
                    className={`rounded-full border px-3 py-1 text-[10px] uppercase tracking-[0.3em] transition ${
                      overlayEnabled
                        ? 'border-amber-200/60 text-amber-100/90'
                        : 'border-white/10 text-white/50 hover:text-white/70'
                    }`}
                    onClick={() => setOverlayEnabled((prev) => !prev)}
                  >
                    Overlay {overlayEnabled ? 'On' : 'Off'}
                  </button>
                  <button
                    className="rounded-full border border-white/10 bg-white/5 p-1.5 text-white/60 transition hover:text-white"
                    onClick={handleAddStructure}
                  >
                    <Plus className="h-4 w-4" />
                  </button>
                </div>
              </div>
              <div className="mt-4 space-y-3">
                {structures.map((structure) => (
                  <div
                    key={structure.id}
                    className={`group cursor-pointer rounded-xl border bg-slate-900/60 p-3 transition ${
                      structure.id === activeId
                        ? 'border-amber-200/60 shadow-lg shadow-amber-400/20'
                        : 'border-white/10 hover:border-white/30'
                    }`}
                    onClick={() => handleSelectStructure(structure.id)}
                  >
                    <div className="flex items-center justify-between">
                      <div className="flex items-center gap-3">
                        <div
                          className={`h-10 w-10 rounded-xl bg-gradient-to-br ${structure.color}`}
                        />
                        <div>
                          <p className="text-sm font-semibold">{structure.name}</p>
                          <p className="text-xs text-white/50">ID {structure.id}</p>
                        </div>
                      </div>
                      <button className="text-white/40 transition group-hover:text-white">
                        <ChevronDown className="h-4 w-4" />
                      </button>
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
                      className="mt-3 space-y-2 text-xs text-white/60"
                      onClick={(event) => event.stopPropagation()}
                    >
                      <label className="flex items-center justify-between gap-2">
                        <span>Visible</span>
                        <input
                          type="checkbox"
                          className="h-4 w-4 accent-amber-300"
                          checked={structure.isVisible}
                          onChange={(event) =>
                            handleToggleVisibility(structure.id, event.target.checked)
                          }
                        />
                      </label>
                      <div className="flex items-center justify-between gap-2">
                        <span>Opacity</span>
                        <div className="flex items-center gap-2">
                          <input
                            type="range"
                            min={0}
                            max={100}
                            step={1}
                            value={Math.round(structure.opacity * 100)}
                            onChange={(event) =>
                              handleOpacityChange(
                                structure.id,
                                Number(event.target.value) / 100,
                              )
                            }
                            className="h-1 w-24 cursor-pointer accent-amber-300"
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
              <div className="mt-4 grid gap-2">
                <button
                  className="flex w-full items-center justify-center gap-2 rounded-xl border border-white/10 bg-white/5 px-3 py-2 text-xs text-white/70 transition hover:border-white/30"
                  onClick={() => fileInputRef.current?.click()}
                >
                  <ArrowDownToLine className="h-4 w-4" />
                  Import .in File
                </button>
                <button
                  className="flex w-full items-center justify-center gap-2 rounded-xl border border-white/10 bg-white/5 px-3 py-2 text-xs text-white/70 transition hover:border-white/30"
                  onClick={handleImportClipboard}
                >
                  Import from Clipboard
                </button>
                <input
                  ref={fileInputRef}
                  type="file"
                  accept=".in,.txt"
                  className="hidden"
                  onChange={handleImportFile}
                />
              </div>
            </div>
            <div className="rounded-2xl border border-white/10 bg-white/5 p-4">
              <div className="flex items-center justify-between">
                <p className="text-xs uppercase tracking-[0.3em] text-white/50">
                  Compare
                </p>
                <RefreshCw className="h-4 w-4 text-white/40" />
              </div>
              <div className="mt-4 space-y-3 text-sm text-white/70">
                <div className="flex items-center justify-between gap-3">
                  <span>Compare With</span>
                  <select
                    value={compareTarget?.id ?? ''}
                    onChange={(event) => setCompareTargetId(event.target.value || null)}
                    className="rounded-lg border border-white/10 bg-slate-950/70 px-2 py-1 text-xs text-white/80 focus:border-amber-300 focus:outline-none"
                  >
                    {compareCandidates.length === 0 ? (
                      <option value="">None</option>
                    ) : null}
                    {compareCandidates.map((candidate) => (
                      <option key={candidate.id} value={candidate.id}>
                        {candidate.name}
                      </option>
                    ))}
                  </select>
                </div>
                {countMismatch ? (
                  <div className="rounded-lg border border-rose-400/30 bg-rose-500/10 px-3 py-2 text-xs text-rose-200">
                    原子数が一致しません。インデックス対応の転写・距離に影響します。
                  </div>
                ) : null}
                <div className="flex items-center justify-between">
                  <span>PBC</span>
                  <span className="rounded-full border border-white/10 px-3 py-1 text-xs">
                    {pbcState.enabled ? 'On' : 'Off'}
                  </span>
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
                              #{row.index + 1} {row.source.symbol}→
                              {row.target.symbol}
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
                    {selectedIndices.length} / {atoms.length}
                  </span>
                </div>
              </div>
              <button className="mt-4 w-full rounded-xl bg-white/10 px-3 py-2 text-xs text-white transition hover:bg-white/20">
                Align Selected
              </button>
              <button
                className="mt-2 w-full rounded-xl border border-amber-200/40 px-3 py-2 text-xs text-amber-100/90 transition hover:text-amber-200 disabled:cursor-not-allowed disabled:opacity-40"
                onClick={handleSurfaceTransfer}
                disabled={!compareTarget}
              >
                Surface Transfer ({activeId}→{compareTarget?.id ?? '—'})
              </button>
            </div>
          </section>

          <section className="flex w-full flex-1 flex-col gap-4">
            <div className="rounded-2xl border border-white/10 bg-white/5 p-4 shadow-lg shadow-black/40">
              <div className="flex flex-wrap items-center justify-between gap-3">
                <div>
                  <p className="text-sm uppercase tracking-[0.3em] text-white/50">
                    Atom Table
                  </p>
                  <p className="text-lg font-semibold">Editable Coordinates</p>
                </div>
                <div className="flex flex-wrap items-center gap-2 text-xs text-white/60">
                  <span className="rounded-full border border-white/10 px-3 py-1">
                    XYZ mode
                  </span>
                  <span className="rounded-full border border-white/10 px-3 py-1">
                    Ångström
                  </span>
                  <button
                    className="rounded-full border border-white/10 px-3 py-1 text-white/70 transition hover:text-white"
                    onClick={handleExport}
                  >
                    Export .in
                  </button>
                  <button
                    className="rounded-full border border-amber-200/40 px-3 py-1 text-amber-100/90 transition hover:text-amber-200"
                    onClick={handleShareHtml}
                  >
                    Share HTML
                  </button>
                </div>
              </div>
              <div className="mt-4 rounded-xl border border-white/10 bg-slate-900/60 p-4">
                <div className="flex flex-wrap items-center justify-between gap-2 text-xs text-white/60">
                  <span className="uppercase tracking-[0.3em]">Input</span>
                  <div className="flex items-center gap-2">
                    <button
                      className="rounded-full border border-white/10 px-3 py-1 text-white/70 transition hover:text-white"
                      onClick={() => setQeInput(sampleQe)}
                    >
                      Reset
                    </button>
                    <button
                      className="rounded-full bg-amber-300 px-4 py-1 text-slate-900 transition hover:bg-amber-200"
                      onClick={handleParse}
                      disabled={isParsing}
                    >
                      {isParsing ? 'Parsing...' : 'Parse'}
                    </button>
                  </div>
                </div>
                <textarea
                  value={qeInput}
                  onChange={(event) => setQeInput(event.target.value)}
                  className="mt-3 h-32 w-full resize-none rounded-lg border border-white/10 bg-slate-950/80 p-3 font-mono text-xs text-white/80 focus:border-amber-300 focus:outline-none"
                  aria-label="Quantum ESPRESSO input"
                />
                {error ? (
                  <p className="mt-2 text-xs text-rose-300">{error}</p>
                ) : null}
                {exported ? (
                  <p className="mt-2 text-xs text-emerald-200">
                    エクスポート済み（クリップボードにコピーしました）
                  </p>
                ) : null}
              </div>
              <div className="mt-4 overflow-hidden rounded-xl border border-white/10">
                <div className="flex items-center justify-between border-b border-white/10 bg-slate-900/70 px-3 py-2 text-xs text-white/60">
                  <span>Atoms</span>
                  <button
                    className="rounded-full border border-white/10 px-3 py-1 text-xs text-white/70 transition hover:text-white"
                    onClick={handleAddAtom}
                  >
                    Add Atom
                  </button>
                </div>
                <div className="max-h-[360px] overflow-auto bg-slate-950/60">
                  <table className="w-full text-left text-xs">
                    <thead className="sticky top-0 bg-slate-950/90 text-white/60">
                      <tr>
                        {['Sel', '#', 'Atom', 'x', 'y', 'z'].map((label) => (
                          <th key={label} className="px-3 py-2 font-medium">
                            {label}
                          </th>
                        ))}
                      </tr>
                    </thead>
                    <tbody>
                      {atoms.length === 0 ? (
                        <tr>
                          <td
                            colSpan={6}
                            className="px-3 py-6 text-center text-white/50"
                          >
                            まだ構造が読み込まれていません
                          </td>
                        </tr>
                      ) : null}
                      {atoms.map((atom, index) => {
                        const draft =
                          atomDrafts[index] ?? ({
                            symbol: atom.symbol,
                            x: atom.x.toFixed(4),
                            y: atom.y.toFixed(4),
                            z: atom.z.toFixed(4),
                          } satisfies AtomDraft)
                        return (
                          <Fragment key={`${atom.symbol}-${index}`}>
                            <tr className="border-t border-white/5">
                              <td className="px-3 py-2">
                                <input
                                  type="checkbox"
                                  checked={selectedIndices.includes(index)}
                                  onChange={() => toggleSelect(index)}
                                  className="h-4 w-4 rounded border-white/20 bg-slate-900/70 text-amber-300 accent-amber-300"
                                  aria-label={`select atom ${index + 1}`}
                                />
                              </td>
                              <td className="px-3 py-2 text-white/40">
                                {index + 1}
                              </td>
                              <td className="px-3 py-2">
                                <input
                                  value={draft.symbol}
                                  onChange={(event) =>
                                    handleAtomChange(
                                      index,
                                      'symbol',
                                      event.target.value,
                                    )
                                  }
                                  className="w-16 rounded-md border border-white/10 bg-slate-900/70 px-2 py-1 text-white/80 focus:border-amber-300 focus:outline-none"
                                  aria-label={`atom ${index + 1} symbol`}
                                />
                              </td>
                              {(['x', 'y', 'z'] as const).map((axis) => (
                                <td key={axis} className="px-3 py-2">
                                  <input
                                    value={draft[axis]}
                                    onChange={(event) =>
                                      handleAtomChange(
                                        index,
                                        axis,
                                        event.target.value,
                                      )
                                    }
                                    className="w-24 rounded-md border border-white/10 bg-slate-900/70 px-2 py-1 text-white/80 focus:border-amber-300 focus:outline-none"
                                    aria-label={`atom ${index + 1} ${axis}`}
                                  />
                                </td>
                              ))}
                            </tr>
                          </Fragment>
                        )
                      })}
                    </tbody>
                  </table>
                </div>
              </div>
            </div>

            <div className="grid gap-4 lg:grid-cols-2">
              <div className="rounded-2xl border border-white/10 bg-white/5 p-4">
                <div className="flex items-center justify-between">
                  <p className="text-xs uppercase tracking-[0.3em] text-white/50">
                    Selection Tools
                  </p>
                  <PencilRuler className="h-4 w-4 text-white/40" />
                </div>
                <div className="mt-4 grid gap-3 text-sm">
                  <div className="rounded-xl border border-white/10 bg-slate-900/60 p-3">
                    <p className="text-white/60">Shift by vector</p>
                    <div className="mt-2 grid grid-cols-3 gap-2 text-xs">
                      {(['x', 'y', 'z'] as const).map((axis) => (
                        <input
                          key={axis}
                          value={shiftDraft[axis]}
                          onChange={(event) =>
                            setShiftDraft((prev) => ({
                              ...prev,
                              [axis]: event.target.value,
                            }))
                          }
                          className="rounded-lg border border-white/10 bg-white/5 px-2 py-2 text-white/70 focus:border-amber-300 focus:outline-none"
                          placeholder={`d${axis}`}
                        />
                      ))}
                    </div>
                  </div>
                  <div className="grid grid-cols-2 gap-2 text-xs">
                    <button
                      className="rounded-xl bg-white/10 px-3 py-2 text-white transition hover:bg-white/20"
                      onClick={handleShiftSelected}
                    >
                      Apply Shift
                    </button>
                    <button
                      className="rounded-xl border border-white/10 px-3 py-2 text-white/70 transition hover:text-white"
                      onClick={handleAlignOrigin}
                    >
                      Align to Origin
                    </button>
                    <button
                      className="rounded-xl border border-white/10 px-3 py-2 text-white/70 transition hover:text-white"
                      onClick={handleAlignCentroid}
                    >
                      Align to Centroid
                    </button>
                    <button
                      className="rounded-xl border border-white/10 px-3 py-2 text-white/70 transition hover:text-white"
                      onClick={handleCopySelected}
                    >
                      Copy Selected
                    </button>
                    <button
                      className="rounded-xl border border-white/10 px-3 py-2 text-white/70 transition hover:text-white"
                      onClick={handlePasteAppend}
                    >
                      Paste Append
                    </button>
                    <button
                      className="rounded-xl border border-white/10 px-3 py-2 text-white/70 transition hover:text-white"
                      onClick={() => setSelectedIndices([])}
                    >
                      Clear Selection
                    </button>
                  </div>
                </div>
              </div>
              <div className="rounded-2xl border border-white/10 bg-white/5 p-4">
                <div className="flex items-center justify-between">
                  <p className="text-xs uppercase tracking-[0.3em] text-white/50">
                    View Controls
                  </p>
                  <Droplets className="h-4 w-4 text-white/40" />
                </div>
                <div className="mt-4 space-y-3 text-sm text-white/70">
                  <div className="flex items-center justify-between">
                    <span>Overlay Mode</span>
                    <span className="rounded-full border border-white/10 px-3 py-1 text-xs">
                      {overlayEnabled ? 'On' : 'Off'}
                    </span>
                  </div>
                  <div className="flex items-center justify-between">
                    <span>Visible Structures</span>
                    <span className="rounded-full border border-white/10 px-3 py-1 text-xs">
                      {visibleOverlayCount}
                    </span>
                  </div>
                  <div className="flex items-center justify-between">
                    <span>Active Opacity</span>
                    <span className="text-xs">
                      {Math.round((activeStructure?.opacity ?? 1) * 100)}%
                    </span>
                  </div>
                </div>
              </div>
            </div>
          </section>

          <section className="flex w-full flex-col gap-4 lg:max-w-[360px]">
            <div className="flex flex-1 flex-col rounded-2xl border border-white/10 bg-white/5 p-4 shadow-lg shadow-black/40">
              <div className="flex items-center justify-between">
                <div>
                  <p className="text-xs uppercase tracking-[0.3em] text-white/50">
                    3D Viewer
                  </p>
                  <p className="text-lg font-semibold">Mol* Preview</p>
                </div>
                <Layers className="h-4 w-4 text-white/40" />
              </div>
              <div className="mt-4 h-80 flex-1 lg:h-full">
                <MolstarViewer structures={viewerStructures} />
              </div>
              <div className="mt-4 grid grid-cols-2 gap-2 text-xs text-white/70">
                <button className="rounded-lg border border-white/10 px-3 py-2 transition hover:border-white/30">
                  Auto Fit
                </button>
                <button className="rounded-lg border border-white/10 px-3 py-2 transition hover:border-white/30">
                  Reset View
                </button>
                <button className="rounded-lg border border-white/10 px-3 py-2 transition hover:border-white/30">
                  Snapshot
                </button>
                <button className="rounded-lg border border-white/10 px-3 py-2 transition hover:border-white/30">
                  Clip
                </button>
              </div>
            </div>
          </section>
        </main>
      </div>
    </div>
  )
}
