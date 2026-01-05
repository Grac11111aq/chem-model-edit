import { createFileRoute } from '@tanstack/react-router'
import { useMemo, useState } from 'react'
import { Boxes, Grid3X3, Layers, Wand2 } from 'lucide-react'
import MolstarViewer from '../components/molstar/MolstarViewer'
import { generateSupercell } from '../lib/api'
import { atomsToPdb } from '../lib/pdb'
import type { Atom } from '../lib/types'
import { parseXyzBlock } from '../lib/xyz'

export const Route = createFileRoute('/supercell')({
  component: SupercellPage,
})

function SupercellPage() {
  const sampleXyzA = ['O 0 0 0', 'H 0.757 0.586 0', 'H -0.757 0.586 0'].join(
    '\n',
  )
  const sampleXyzB = ['C 0 0 0', 'H 0.629 0.629 0.629'].join('\n')
  const [xyzA, setXyzA] = useState(sampleXyzA)
  const [xyzB, setXyzB] = useState(sampleXyzB)
  const [sequence, setSequence] = useState('AB,AA')
  const [lattice, setLattice] = useState({
    a: { x: '5.0', y: '0.0', z: '0.0' },
    b: { x: '0.0', y: '5.0', z: '0.0' },
    c: { x: '0.0', y: '0.0', z: '5.0' },
  })
  const [atoms, setAtoms] = useState<Atom[]>([])
  const [meta, setMeta] = useState({ na: 0, nb: 0, layers: 0 })
  const [error, setError] = useState<string | null>(null)
  const [isGenerating, setIsGenerating] = useState(false)

  const pdb = useMemo(() => atomsToPdb(atoms), [atoms])

  const handleGenerate = async () => {
    setIsGenerating(true)
    setError(null)
    try {
      const structureA = { atoms: parseXyzBlock(xyzA) }
      const structureB = { atoms: parseXyzBlock(xyzB) }
      const parseVector = (vec: { x: string; y: string; z: string }) => ({
        x: Number(vec.x),
        y: Number(vec.y),
        z: Number(vec.z),
      })
      const a = parseVector(lattice.a)
      const b = parseVector(lattice.b)
      const c = parseVector(lattice.c)
      if (
        [a.x, a.y, a.z, b.x, b.y, b.z, c.x, c.y, c.z].some((v) =>
          Number.isNaN(v),
        )
      ) {
        throw new Error('格子ベクトルの数値が無効です。')
      }
      if (structureA.atoms.length === 0 || structureB.atoms.length === 0) {
        throw new Error('A/B両方の原子が必要です。')
      }
      const result = await generateSupercell({
        structureA,
        structureB,
        sequence,
        lattice: { a, b, c },
      })
      setAtoms(result.structure.atoms)
      setMeta(result.meta)
    } catch (err) {
      setError(err instanceof Error ? err.message : '生成に失敗しました。')
    } finally {
      setIsGenerating(false)
    }
  }

  return (
    <div className="min-h-screen bg-slate-950 text-white">
      <div className="mx-auto max-w-[1200px] px-6 py-10">
        <div className="flex flex-wrap items-center justify-between gap-4">
          <div>
            <p className="text-xs uppercase tracking-[0.3em] text-white/50">
              Supercell Builder
            </p>
            <h1 className="mt-2 text-3xl font-semibold">Stack slabs with intent</h1>
            <p className="mt-2 text-white/60">
              Define lattice parameters and stacking sequences to generate a
              supercell preview.
            </p>
          </div>
          <button
            className="flex items-center gap-2 rounded-full bg-amber-300 px-4 py-2 text-sm font-semibold text-slate-900 shadow-lg shadow-amber-400/30"
            onClick={handleGenerate}
            disabled={isGenerating}
          >
            <Wand2 className="h-4 w-4" />
            {isGenerating ? 'Generating...' : 'Generate'}
          </button>
        </div>

        <div className="mt-8 grid gap-6 lg:grid-cols-[1.2fr_0.8fr]">
          <div className="rounded-2xl border border-white/10 bg-white/5 p-6">
            <div className="flex items-center gap-3 text-white/70">
              <Grid3X3 className="h-5 w-5" />
              <span className="text-sm uppercase tracking-[0.3em]">Lattice</span>
            </div>
            <div className="mt-6 grid grid-cols-3 gap-3 text-sm">
              {(['a', 'b', 'c'] as const).map((axis) => (
                <div
                  key={axis}
                  className="rounded-xl border border-white/10 bg-slate-900/60 px-3 py-4"
                >
                  <p className="text-xs text-white/40">{axis}-vector</p>
                  <div className="mt-2 space-y-1 text-xs text-white/60">
                    {(['x', 'y', 'z'] as const).map((component) => (
                      <input
                        key={`${axis}-${component}`}
                        value={lattice[axis][component]}
                        onChange={(event) =>
                          setLattice((prev) => ({
                            ...prev,
                            [axis]: {
                              ...prev[axis],
                              [component]: event.target.value,
                            },
                          }))
                        }
                        className="w-full rounded-md border border-white/10 bg-slate-950/70 px-2 py-1 text-white/80 focus:border-amber-300 focus:outline-none"
                      />
                    ))}
                  </div>
                </div>
              ))}
            </div>
            <div className="mt-6 rounded-xl border border-white/10 bg-slate-900/60 p-4">
              <p className="text-xs uppercase tracking-[0.3em] text-white/50">
                Sequence
              </p>
              <input
                value={sequence}
                onChange={(event) => setSequence(event.target.value)}
                className="mt-3 w-full rounded-md border border-white/10 bg-slate-950/70 px-3 py-2 text-sm text-white/80 focus:border-amber-300 focus:outline-none"
              />
            </div>
            <div className="mt-6 grid gap-4 lg:grid-cols-2">
              <div className="rounded-xl border border-white/10 bg-slate-900/60 p-4">
                <p className="text-xs uppercase tracking-[0.3em] text-white/50">
                  Structure A (XYZ)
                </p>
                <textarea
                  value={xyzA}
                  onChange={(event) => setXyzA(event.target.value)}
                  className="mt-3 h-28 w-full resize-none rounded-md border border-white/10 bg-slate-950/70 p-2 font-mono text-xs text-white/80 focus:border-amber-300 focus:outline-none"
                  aria-label="Structure A XYZ"
                />
              </div>
              <div className="rounded-xl border border-white/10 bg-slate-900/60 p-4">
                <p className="text-xs uppercase tracking-[0.3em] text-white/50">
                  Structure B (XYZ)
                </p>
                <textarea
                  value={xyzB}
                  onChange={(event) => setXyzB(event.target.value)}
                  className="mt-3 h-28 w-full resize-none rounded-md border border-white/10 bg-slate-950/70 p-2 font-mono text-xs text-white/80 focus:border-amber-300 focus:outline-none"
                  aria-label="Structure B XYZ"
                />
              </div>
            </div>
            {error ? <p className="mt-4 text-xs text-rose-300">{error}</p> : null}
          </div>

          <div className="flex flex-col gap-6">
            <div className="rounded-2xl border border-white/10 bg-white/5 p-6">
              <div className="flex items-center gap-3 text-white/70">
                <Boxes className="h-5 w-5" />
                <span className="text-sm uppercase tracking-[0.3em]">Output</span>
              </div>
              <div className="mt-6 space-y-3 text-sm text-white/70">
                <div className="flex items-center justify-between">
                  <span>Supercell size</span>
                  <span className="font-semibold text-white">
                    {meta.na} × {meta.nb} × 1
                  </span>
                </div>
                <div className="flex items-center justify-between">
                  <span>Atoms</span>
                  <span className="font-semibold text-white">
                    {atoms.length}
                  </span>
                </div>
                <div className="flex items-center justify-between">
                  <span>Layers</span>
                  <span className="font-semibold text-white">{meta.layers}</span>
                </div>
              </div>
            </div>
            <div className="flex flex-1 flex-col rounded-2xl border border-white/10 bg-white/5 p-6">
              <div className="flex items-center justify-between text-white/60">
                <span className="text-xs uppercase tracking-[0.3em]">Viewer</span>
                <Layers className="h-4 w-4" />
              </div>
              <div className="mt-4 min-h-[240px] flex-1">
                <MolstarViewer pdbText={atoms.length ? pdb : undefined} />
              </div>
            </div>
          </div>
        </div>
      </div>
    </div>
  )
}
