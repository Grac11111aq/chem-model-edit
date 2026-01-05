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
import type { Atom } from '../lib/types'
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
    },
    {
      id: 'B',
      name: 'Ligand-B',
      color: palette[1],
      atoms: [],
      drafts: [],
    },
  ])
  const [activeId, setActiveId] = useState(structures[0]?.id ?? 'A')
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
  const displayPdb = useMemo(() => {
    if (atoms.length === 0) {
      return samplePdb
    }
    return atomsToPdb(atoms)
  }, [atoms, samplePdb])
  const selectedAtoms = useMemo(
    () => selectedIndices.map((index) => atoms[index]).filter(Boolean),
    [atoms, selectedIndices],
  )
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

  const syncDrafts = (nextAtoms: Atom[]) => {
    updateActive((structure) => ({
      ...structure,
      drafts: nextAtoms.map((atom) => ({
        symbol: atom.symbol,
        x: atom.x.toFixed(4),
        y: atom.y.toFixed(4),
        z: atom.z.toFixed(4),
      })),
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
        drafts: structure.atoms.map((atom) => ({
          symbol: atom.symbol,
          x: atom.x.toFixed(4),
          y: atom.y.toFixed(4),
          z: atom.z.toFixed(4),
        })),
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
    if (atoms.length === 0) {
      setError('共有する原子がありません。')
      return
    }
    setError(null)
    downloadShareHtml(atoms)
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
      drafts: [
        ...structure.drafts,
        { symbol: 'X', x: '0.0000', y: '0.0000', z: '0.0000' },
      ],
    }))
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
        drafts: [
          ...structure.drafts,
          ...parsed.map((atom) => ({
            symbol: atom.symbol,
            x: atom.x.toFixed(4),
            y: atom.y.toFixed(4),
            z: atom.z.toFixed(4),
          })),
        ],
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
                <button
                  className="rounded-full border border-white/10 bg-white/5 p-1.5 text-white/60 transition hover:text-white"
                  onClick={handleAddStructure}
                >
                  <Plus className="h-4 w-4" />
                </button>
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
                <div className="flex items-center justify-between">
                  <span>RMSD</span>
                  <span className="font-semibold text-white">0.042 Å</span>
                </div>
                <div className="flex items-center justify-between">
                  <span>Pair Distance</span>
                  <span className="font-semibold text-white">
                    {selectedDistance ? `${selectedDistance.toFixed(4)} Å` : '—'}
                  </span>
                </div>
                <div className="rounded-xl border border-white/10 bg-slate-900/60 px-3 py-2 text-xs text-white/60">
                  <p className="uppercase tracking-[0.3em] text-white/40">
                    Selected Distances
                  </p>
                  <div className="mt-2 space-y-1">
                    {selectedAtoms.length < 2 ? (
                      <p className="text-white/40">2つ以上選択してください。</p>
                    ) : (
                      selectedAtoms.slice(1).map((atom, index) => {
                        const target = selectedAtoms[0]
                        const dx = atom.x - target.x
                        const dy = atom.y - target.y
                        const dz = atom.z - target.z
                        const dist = Math.sqrt(dx * dx + dy * dy + dz * dz)
                        return (
                          <div key={`${atom.symbol}-${index}`} className="flex justify-between">
                            <span>
                              {target.symbol}→{atom.symbol}
                            </span>
                            <span className="text-white">{dist.toFixed(4)} Å</span>
                          </div>
                        )
                      })
                    )}
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
                    <span>Background</span>
                    <span className="rounded-full border border-white/10 px-3 py-1 text-xs">
                      Polar Dark
                    </span>
                  </div>
                  <div className="flex items-center justify-between">
                    <span>Labels</span>
                    <span className="rounded-full border border-amber-200/50 px-3 py-1 text-xs text-amber-200">
                      On
                    </span>
                  </div>
                  <div className="flex items-center justify-between">
                    <span>Opacity</span>
                    <span className="text-xs">65%</span>
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
                <MolstarViewer pdbText={displayPdb} />
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
