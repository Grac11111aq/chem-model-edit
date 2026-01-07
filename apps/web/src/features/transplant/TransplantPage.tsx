import { useRef, useState } from 'react'
import { Clipboard, FileText, Play, Upload } from 'lucide-react'

import { Button } from '@/components/ui/button'
import { Card } from '@/components/ui/card'
import { Textarea } from '@/components/ui/textarea'
import { deltaTransplant } from '@/lib/api'

export default function TransplantPage() {
  const [inputs, setInputs] = useState({
    smallIn: '',
    smallOut: '',
    largeIn: '',
  })
  const [result, setResult] = useState('')
  const [error, setError] = useState<string | null>(null)
  const [isRunning, setIsRunning] = useState(false)

  const smallInRef = useRef<HTMLInputElement | null>(null)
  const smallOutRef = useRef<HTMLInputElement | null>(null)
  const largeInRef = useRef<HTMLInputElement | null>(null)

  const updateInput = (key: keyof typeof inputs, value: string) => {
    setInputs((prev) => ({ ...prev, [key]: value }))
  }

  const handleFile = async (
    event: React.ChangeEvent<HTMLInputElement>,
    key: keyof typeof inputs,
  ) => {
    const file = event.target.files?.[0]
    if (!file) {
      return
    }
    const text = await file.text()
    updateInput(key, text)
    event.target.value = ''
  }

  const handlePaste = async (key: keyof typeof inputs) => {
    try {
      const text = await navigator.clipboard.readText()
      if (!text.trim()) {
        setError('クリップボードが空です。')
        return
      }
      setError(null)
      updateInput(key, text)
    } catch (err) {
      setError(
        err instanceof Error
          ? err.message
          : 'クリップボード読取に失敗しました。',
      )
    }
  }

  const handleRun = async () => {
    if (
      !inputs.smallIn.trim() ||
      !inputs.smallOut.trim() ||
      !inputs.largeIn.trim()
    ) {
      setError('small .in / small .out / large .in の3つを入力してください。')
      return
    }
    setError(null)
    setIsRunning(true)
    try {
      const content = await deltaTransplant({
        smallIn: inputs.smallIn,
        smallOut: inputs.smallOut,
        largeIn: inputs.largeIn,
      })
      setResult(content)
    } catch (err) {
      setResult('')
      setError(err instanceof Error ? err.message : 'Δ移植に失敗しました。')
    } finally {
      setIsRunning(false)
    }
  }

  const handleCopyResult = async () => {
    if (!result) {
      setError('コピーする結果がありません。')
      return
    }
    setError(null)
    try {
      await navigator.clipboard.writeText(result)
    } catch (err) {
      setError(err instanceof Error ? err.message : 'コピーに失敗しました。')
    }
  }

  const handleDownloadResult = () => {
    if (!result) {
      setError('ダウンロードする結果がありません。')
      return
    }
    setError(null)
    const blob = new Blob([result], { type: 'text/plain' })
    const url = URL.createObjectURL(blob)
    const a = document.createElement('a')
    a.href = url
    a.download = 'transplanted.in'
    a.click()
    URL.revokeObjectURL(url)
  }

  return (
    <div className="min-h-screen bg-slate-950 text-white">
      <div className="relative overflow-hidden">
        <div className="pointer-events-none absolute inset-0 bg-[radial-gradient(circle_at_top,_rgba(94,234,212,0.18),_transparent_55%)]" />
        <div className="pointer-events-none absolute inset-0 bg-[radial-gradient(circle_at_15%_30%,_rgba(248,113,113,0.12),_transparent_45%)]" />
        <div className="pointer-events-none absolute inset-0 bg-[radial-gradient(circle_at_85%_20%,_rgba(59,130,246,0.16),_transparent_45%)]" />
        <main className="relative mx-auto flex min-h-[calc(100vh-64px)] max-w-[1400px] flex-col gap-6 px-6 py-6 lg:flex-row">
          <section className="flex w-full flex-1 flex-col gap-4">
            <Card className="border-white/10 bg-white/5 p-5 shadow-lg shadow-black/40">
              <div className="flex flex-wrap items-center justify-between gap-3">
                <div>
                  <p className="text-xs uppercase tracking-[0.3em] text-white/50">
                    Delta Transplant
                  </p>
                  <p className="text-lg font-semibold">
                    小スラブの緩和変位を大スラブへ
                  </p>
                </div>
                <Button
                  className="rounded-full bg-emerald-300 text-xs font-semibold text-slate-900 hover:bg-emerald-200"
                  onClick={handleRun}
                  disabled={isRunning}
                >
                  <Play className="h-4 w-4" />
                  {isRunning ? 'Running...' : 'Run Transplant'}
                </Button>
              </div>
              <div className="mt-4 rounded-xl border border-white/10 bg-slate-900/60 px-4 py-3 text-xs text-white/60">
                <p className="uppercase tracking-[0.3em] text-white/40">
                  Notes
                </p>
                <ul className="mt-2 space-y-1">
                  <li>ATOMIC_POSITIONS の 0/1 フラグ必須</li>
                  <li>非直交セルは未対応（wrap は Lx/Ly のみ）</li>
                  <li>可動原子数が一致しない場合はエラー</li>
                </ul>
              </div>
            </Card>

            <div className="grid gap-4 lg:grid-cols-3">
              <Card className="border-white/10 bg-white/5 p-4">
                <div className="flex items-center justify-between text-xs text-white/60">
                  <span className="uppercase tracking-[0.3em]">Small .in</span>
                  <div className="flex gap-2">
                    <Button
                      variant="outline"
                      size="sm"
                      className="rounded-full border-white/10 bg-white/5 text-[10px] uppercase tracking-[0.2em] text-white/70 hover:text-white"
                      onClick={() => smallInRef.current?.click()}
                    >
                      <Upload className="inline h-3 w-3" /> File
                    </Button>
                    <Button
                      variant="outline"
                      size="sm"
                      className="rounded-full border-white/10 bg-white/5 text-[10px] uppercase tracking-[0.2em] text-white/70 hover:text-white"
                      onClick={() => handlePaste('smallIn')}
                    >
                      <Clipboard className="inline h-3 w-3" /> Paste
                    </Button>
                  </div>
                </div>
                <Textarea
                  value={inputs.smallIn}
                  onChange={(event) =>
                    updateInput('smallIn', event.target.value)
                  }
                  className="mt-3 h-64 resize-none border-white/10 bg-slate-950/80 font-mono text-xs text-white/80 focus-visible:ring-emerald-300"
                  placeholder="小スラブ input (.in)"
                />
                <input
                  ref={smallInRef}
                  type="file"
                  accept=".in,.txt"
                  className="hidden"
                  onChange={(event) => handleFile(event, 'smallIn')}
                />
              </Card>

              <Card className="border-white/10 bg-white/5 p-4">
                <div className="flex items-center justify-between text-xs text-white/60">
                  <span className="uppercase tracking-[0.3em]">Small .out</span>
                  <div className="flex gap-2">
                    <Button
                      variant="outline"
                      size="sm"
                      className="rounded-full border-white/10 bg-white/5 text-[10px] uppercase tracking-[0.2em] text-white/70 hover:text-white"
                      onClick={() => smallOutRef.current?.click()}
                    >
                      <Upload className="inline h-3 w-3" /> File
                    </Button>
                    <Button
                      variant="outline"
                      size="sm"
                      className="rounded-full border-white/10 bg-white/5 text-[10px] uppercase tracking-[0.2em] text-white/70 hover:text-white"
                      onClick={() => handlePaste('smallOut')}
                    >
                      <Clipboard className="inline h-3 w-3" /> Paste
                    </Button>
                  </div>
                </div>
                <Textarea
                  value={inputs.smallOut}
                  onChange={(event) =>
                    updateInput('smallOut', event.target.value)
                  }
                  className="mt-3 h-64 resize-none border-white/10 bg-slate-950/80 font-mono text-xs text-white/80 focus-visible:ring-emerald-300"
                  placeholder="小スラブ output (.out)"
                />
                <input
                  ref={smallOutRef}
                  type="file"
                  accept=".out,.txt"
                  className="hidden"
                  onChange={(event) => handleFile(event, 'smallOut')}
                />
              </Card>

              <Card className="border-white/10 bg-white/5 p-4">
                <div className="flex items-center justify-between text-xs text-white/60">
                  <span className="uppercase tracking-[0.3em]">Large .in</span>
                  <div className="flex gap-2">
                    <Button
                      variant="outline"
                      size="sm"
                      className="rounded-full border-white/10 bg-white/5 text-[10px] uppercase tracking-[0.2em] text-white/70 hover:text-white"
                      onClick={() => largeInRef.current?.click()}
                    >
                      <Upload className="inline h-3 w-3" /> File
                    </Button>
                    <Button
                      variant="outline"
                      size="sm"
                      className="rounded-full border-white/10 bg-white/5 text-[10px] uppercase tracking-[0.2em] text-white/70 hover:text-white"
                      onClick={() => handlePaste('largeIn')}
                    >
                      <Clipboard className="inline h-3 w-3" /> Paste
                    </Button>
                  </div>
                </div>
                <Textarea
                  value={inputs.largeIn}
                  onChange={(event) =>
                    updateInput('largeIn', event.target.value)
                  }
                  className="mt-3 h-64 resize-none border-white/10 bg-slate-950/80 font-mono text-xs text-white/80 focus-visible:ring-emerald-300"
                  placeholder="大スラブ input (.in)"
                />
                <input
                  ref={largeInRef}
                  type="file"
                  accept=".in,.txt"
                  className="hidden"
                  onChange={(event) => handleFile(event, 'largeIn')}
                />
              </Card>
            </div>
          </section>

          <section className="flex w-full max-w-[520px] flex-col gap-4">
            <Card className="border-white/10 bg-white/5 p-4 shadow-lg shadow-black/40">
              <div className="flex items-center justify-between">
                <p className="text-xs uppercase tracking-[0.3em] text-white/50">
                  Result
                </p>
                <div className="flex gap-2">
                  <Button
                    variant="outline"
                    size="sm"
                    className="rounded-full border-white/10 bg-white/5 text-[10px] uppercase tracking-[0.2em] text-white/70 hover:text-white"
                    onClick={handleCopyResult}
                  >
                    <Clipboard className="inline h-3 w-3" /> Copy
                  </Button>
                  <Button
                    variant="outline"
                    size="sm"
                    className="rounded-full border-white/10 bg-white/5 text-[10px] uppercase tracking-[0.2em] text-white/70 hover:text-white"
                    onClick={handleDownloadResult}
                  >
                    <FileText className="inline h-3 w-3" /> Download
                  </Button>
                </div>
              </div>
              <Textarea
                value={result}
                readOnly
                className="mt-3 h-[520px] resize-none border-white/10 bg-slate-950/80 font-mono text-xs text-white/80 focus-visible:ring-emerald-300"
                placeholder="移植結果 (.in) がここに表示されます"
              />
              {error ? (
                <p className="mt-2 text-xs text-rose-300">{error}</p>
              ) : null}
              {!error && result ? (
                <p className="mt-2 text-xs text-emerald-200">
                  生成済み（コピー/ダウンロード可）
                </p>
              ) : null}
            </Card>
          </section>
        </main>
      </div>
    </div>
  )
}
