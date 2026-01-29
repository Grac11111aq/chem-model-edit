import { useCallback, useEffect, useMemo, useRef, useState } from 'react'
import {
  Activity,
  AlertTriangle,
  CheckCircle2,
  Copy,
  Download,
  Layers,
  MousePointerClick,
  Play,
  RefreshCw,
  X,
} from 'lucide-react'
import { AtomTable } from './AtomTable'
import { CollapsibleSection } from './CollapsibleSection'
import { SupercellTool } from './SupercellTool'

import type { ToolMode, WorkspaceFile } from '../types'
import type {
  AuthSession,
  Structure,
  SupercellBuildMeta,
  ZPEJobStatus,
  ZPEParseResponse,
  ZPEQueueTarget,
  ZPEResult,
} from '@/lib/types'

import MolstarViewer from '@/components/molstar/MolstarViewer'
import { Badge } from '@/components/ui/badge'
import { Button } from '@/components/ui/button'
import { Input } from '@/components/ui/input'
import {
  Select,
  SelectContent,
  SelectItem,
  SelectTrigger,
  SelectValue,
} from '@/components/ui/select'
import { Switch } from '@/components/ui/switch'
import {
  createEnrollToken,
  createZpeJob,
  downloadZpeFile,
  exportQeInput,
  fetchQueueTargets,
  fetchZpeResult,
  fetchZpeStatus,
  getStructure,
  loginAccount,
  logoutAccount,
  parseZpeInput,
  registerAccount,
  selectQueueTarget,
  structureViewUrl,
} from '@/lib/api'
import { clearSession, getStoredSession, storeSession } from '@/lib/auth'
import { atomsToPdb } from '@/lib/pdb'
import { cn } from '@/lib/utils'

interface ToolPanelProps {
  mode: ToolMode
  files?: Array<WorkspaceFile>
  onClose?: () => void
  variant?: 'stack' | 'dock'
  showHeader?: boolean
  showClose?: boolean
  className?: string
  structures?: Array<WorkspaceFile>
  onSupercellCreated?: (result: {
    structureId: string
    meta: SupercellBuildMeta
  }) => void
}

const toolTitles: Record<ToolMode, string> = {
  transfer: 'Structure Transfer',
  supercell: 'Supercell Builder',
  vibration: 'ZPE / Vibrations',
}

type TransferSummary = {
  structure: Structure
  pdbText: string
  sourceAtoms: number
  targetAtoms: number
  transferredAtoms: number
}

const formatNumber = (value: number | null | undefined, digits = 3) => {
  if (value === null || value === undefined) {
    return '—'
  }
  return value.toFixed(digits)
}

const statusTone = (status?: string | null) => {
  switch (status) {
    case 'finished':
      return 'border-emerald-200 bg-emerald-50 text-emerald-700'
    case 'failed':
      return 'border-red-200 bg-red-50 text-red-700'
    case 'started':
      return 'border-blue-200 bg-blue-50 text-blue-700'
    case 'queued':
      return 'border-amber-200 bg-amber-50 text-amber-700'
    default:
      return 'border-slate-200 bg-slate-100 text-slate-600'
  }
}

const statusLabel = (status?: string | null) => {
  if (!status) {
    return 'idle'
  }
  return status
}

const downloadTextFile = (content: string, filename: string, type: string) => {
  const blob = new Blob([content], { type })
  const url = URL.createObjectURL(blob)
  const link = document.createElement('a')
  link.href = url
  link.download = filename
  link.click()
  URL.revokeObjectURL(url)
}

const createTransferFilename = (targetName: string) => {
  const trimmed = targetName.trim()
  if (!trimmed) {
    return 'Transferred.in'
  }
  const match = trimmed.match(/\.[^/.]+$/)
  const base = match ? trimmed.slice(0, -match[0].length) : trimmed
  const extension = match && match[0].toLowerCase() === '.in' ? match[0] : '.in'
  return `${base}Transferred${extension}`
}

function ZpeToolPanel({ files = [] }: { files?: Array<WorkspaceFile> }) {
  const availableFiles = useMemo(
    () => files.filter((file) => file.kind === 'in'),
    [files],
  )
  const [selectedFileId, setSelectedFileId] = useState<string | null>(null)
  const [parseResult, setParseResult] = useState<ZPEParseResponse | null>(null)
  const [parseError, setParseError] = useState<string | null>(null)
  const [isParsing, setIsParsing] = useState(false)
  const [viewerError, setViewerError] = useState<string | null>(null)
  const [mobileIndices, setMobileIndices] = useState<Set<number>>(new Set())
  const [atomFilter, setAtomFilter] = useState('')
  const [jobId, setJobId] = useState<string | null>(null)
  const [jobStatus, setJobStatus] = useState<ZPEJobStatus | null>(null)
  const [jobResult, setJobResult] = useState<ZPEResult | null>(null)
  const [runError, setRunError] = useState<string | null>(null)
  const [isSubmitting, setIsSubmitting] = useState(false)
  const [useEnviron, setUseEnviron] = useState(false)
  const [calcMode, setCalcMode] = useState<'new' | 'continue'>('continue')
  const [inputDir, setInputDir] = useState('')
  const [session, setSession] = useState<AuthSession | null>(null)
  const [authMode, setAuthMode] = useState<'login' | 'register'>('login')
  const [authEmail, setAuthEmail] = useState('')
  const [authPassword, setAuthPassword] = useState('')
  const [authError, setAuthError] = useState<string | null>(null)
  const [authBusy, setAuthBusy] = useState(false)
  const [queueTargets, setQueueTargets] = useState<Array<ZPEQueueTarget>>([])
  const [activeTargetId, setActiveTargetId] = useState<string | null>(null)
  const [targetsError, setTargetsError] = useState<string | null>(null)
  const [targetsBusy, setTargetsBusy] = useState(false)
  const [enrollToken, setEnrollToken] = useState<{
    token: string
    expires_at: string
    ttl_seconds: number
    label?: string | null
  } | null>(null)
  const [enrollError, setEnrollError] = useState<string | null>(null)
  const [enrollBusy, setEnrollBusy] = useState(false)
  const [tokenCopied, setTokenCopied] = useState(false)
  const parseTokenRef = useRef(0)
  const sessionTokenRef = useRef<string | null>(null)

  useEffect(() => {
    setSession(getStoredSession())
  }, [])

  useEffect(() => {
    sessionTokenRef.current = session?.token ?? null
  }, [session])

  const refreshTargets = useCallback(async () => {
    if (!session) {
      return
    }
    const tokenSnapshot = session.token
    setTargetsBusy(true)
    setTargetsError(null)
    try {
      const payload = await fetchQueueTargets()
      if (sessionTokenRef.current !== tokenSnapshot) {
        return
      }
      setQueueTargets(payload.targets)
      setActiveTargetId(payload.active_target_id ?? null)
    } catch (err) {
      if (sessionTokenRef.current !== tokenSnapshot) {
        return
      }
      setTargetsError(
        err instanceof Error ? err.message : 'Failed to load compute targets.',
      )
    } finally {
      if (sessionTokenRef.current === tokenSnapshot) {
        setTargetsBusy(false)
      }
    }
  }, [session])

  useEffect(() => {
    if (!session) {
      setQueueTargets([])
      setActiveTargetId(null)
      return
    }
    void refreshTargets()
  }, [refreshTargets, session])

  const activeTarget = useMemo(
    () => queueTargets.find((target) => target.target_id === activeTargetId),
    [queueTargets, activeTargetId],
  )

  useEffect(() => {
    if (availableFiles.length === 0) {
      setSelectedFileId(null)
      return
    }
    if (!selectedFileId) {
      setSelectedFileId(availableFiles[0].id)
      return
    }
    if (!availableFiles.some((file) => file.id === selectedFileId)) {
      setSelectedFileId(availableFiles[0].id)
    }
  }, [availableFiles, selectedFileId])

  const selectedFile = useMemo(
    () =>
      selectedFileId
        ? (availableFiles.find((file) => file.id === selectedFileId) ?? null)
        : null,
    [availableFiles, selectedFileId],
  )

  const baseStructure = useMemo(() => {
    return selectedFile?.structure ?? parseResult?.structure ?? null
  }, [parseResult?.structure, selectedFile?.structure])

  const viewerCifUrl = useMemo(() => {
    if (selectedFile?.cifUrl) {
      return selectedFile.cifUrl
    }
    if (selectedFile?.structureId) {
      return structureViewUrl(selectedFile.structureId, { format: 'cif' })
    }
    return null
  }, [selectedFile?.cifUrl, selectedFile?.structureId])

  useEffect(() => {
    setViewerError(null)
  }, [viewerCifUrl])

  useEffect(() => {
    parseTokenRef.current += 1
    setParseResult(null)
    setParseError(null)
    setViewerError(null)
    setMobileIndices(new Set())
    setJobId(null)
    setJobStatus(null)
    setJobResult(null)
    setRunError(null)
    setAtomFilter('')
  }, [selectedFileId])

  useEffect(() => {
    if (!jobId) {
      return
    }
    const activeRef = { current: true }
    let intervalId: number | null = null

    const isActive = () => activeRef.current

    const pollStatus = async () => {
      try {
        const status = await fetchZpeStatus(jobId)
        if (!isActive()) {
          return
        }
        setJobStatus(status)
        if (status.status === 'finished') {
          const result = await fetchZpeResult(jobId)
          if (!isActive()) {
            return
          }
          setJobResult(result)
          if (intervalId) {
            window.clearInterval(intervalId)
          }
        } else if (status.status === 'failed') {
          if (intervalId) {
            window.clearInterval(intervalId)
          }
        }
      } catch (err) {
        if (!activeRef.current) {
          return
        }
        setRunError(
          err instanceof Error
            ? err.message
            : 'Failed to fetch ZPE job status.',
        )
        if (intervalId) {
          window.clearInterval(intervalId)
        }
      }
    }

    void pollStatus()
    intervalId = window.setInterval(pollStatus, 2000)
    return () => {
      activeRef.current = false
      if (intervalId) {
        window.clearInterval(intervalId)
      }
    }
  }, [jobId])

  const fixedIndexSet = useMemo(
    () => new Set(parseResult?.fixed_indices ?? []),
    [parseResult],
  )

  const mobileIndexList = useMemo(
    () => Array.from(mobileIndices).sort((a, b) => a - b),
    [mobileIndices],
  )

  const disabledAtomIndices = useMemo(
    () => Array.from(fixedIndexSet),
    [fixedIndexSet],
  )

  const atomRows = useMemo(() => {
    if (!baseStructure) {
      return []
    }
    return baseStructure.atoms.map((atom, index) => ({
      index,
      symbol: atom.symbol,
      x: atom.x,
      y: atom.y,
      z: atom.z,
    }))
  }, [baseStructure, fixedIndexSet])

  const filteredRows = useMemo(() => {
    const query = atomFilter.trim().toLowerCase()
    if (!query) {
      return atomRows
    }
    return atomRows.filter((row) => {
      if (row.symbol.toLowerCase().includes(query)) {
        return true
      }
      const indexLabel = String(row.index + 1)
      return indexLabel.includes(query)
    })
  }, [atomFilter, atomRows])

  const speciesEntries = useMemo(
    () => Object.entries(parseResult?.atomic_species ?? {}),
    [parseResult],
  )

  const freqSummary = useMemo(() => {
    if (!jobResult || jobResult.freqs_cm.length === 0) {
      return null
    }
    const values = jobResult.freqs_cm
    const min = Math.min(...values)
    const max = Math.max(...values)
    const imagCount = values.filter((value) => value < 0).length
    return {
      count: values.length,
      min,
      max,
      imagCount,
    }
  }, [jobResult])

  const runParse = useCallback(
    async (content: string, structureId?: string | null) => {
      const token = parseTokenRef.current + 1
      parseTokenRef.current = token
      setIsParsing(true)
      setParseError(null)
      try {
        const result = await parseZpeInput(content, structureId)
        if (parseTokenRef.current !== token) {
          return
        }
        setParseResult(result)
        const fixedSet = new Set(result.fixed_indices)
        const nextMobile = new Set<number>()
        result.structure.atoms.forEach((_atom, index) => {
          if (!fixedSet.has(index)) {
            nextMobile.add(index)
          }
        })
        setMobileIndices(nextMobile)
      } catch (err) {
        if (parseTokenRef.current !== token) {
          return
        }
        setParseError(
          err instanceof Error ? err.message : 'ZPE parse failed to run.',
        )
      } finally {
        if (parseTokenRef.current === token) {
          setIsParsing(false)
        }
      }
    },
    [],
  )

  useEffect(() => {
    if (!selectedFile?.qeInput) {
      return
    }
    void runParse(selectedFile.qeInput, selectedFile.structureId)
  }, [
    runParse,
    selectedFile?.id,
    selectedFile?.qeInput,
    selectedFile?.structureId,
  ])

  const handleAuthSubmit = async () => {
    if (!authEmail.trim() || !authPassword) {
      setAuthError('Email and password are required.')
      return
    }
    setAuthBusy(true)
    setAuthError(null)
    try {
      const payload =
        authMode === 'login'
          ? await loginAccount({
              email: authEmail.trim(),
              password: authPassword,
            })
          : await registerAccount({
              email: authEmail.trim(),
              password: authPassword,
            })
      storeSession(payload)
      setSession(payload)
      setAuthPassword('')
      setEnrollToken(null)
    } catch (err) {
      setAuthError(
        err instanceof Error ? err.message : 'Authentication failed.',
      )
    } finally {
      setAuthBusy(false)
    }
  }

  const handleLogout = async () => {
    setAuthBusy(true)
    try {
      await logoutAccount()
    } catch (_err) {
      // ignore logout errors
    } finally {
      clearSession()
      setSession(null)
      setQueueTargets([])
      setActiveTargetId(null)
      setEnrollToken(null)
      setAuthBusy(false)
    }
  }

  const handleSelectTarget = async (targetId: string) => {
    setTargetsBusy(true)
    setTargetsError(null)
    try {
      const response = await selectQueueTarget(targetId)
      setActiveTargetId(response.active_target_id)
    } catch (err) {
      setTargetsError(
        err instanceof Error ? err.message : 'Failed to select target.',
      )
    } finally {
      setTargetsBusy(false)
    }
  }

  const handleGenerateToken = async () => {
    setEnrollBusy(true)
    setEnrollError(null)
    setTokenCopied(false)
    try {
      const response = await createEnrollToken({ ttlSeconds: 3600 })
      setEnrollToken(response)
    } catch (err) {
      setEnrollError(
        err instanceof Error ? err.message : 'Failed to create enroll token.',
      )
    } finally {
      setEnrollBusy(false)
    }
  }

  const handleCopyToken = async () => {
    if (!enrollToken?.token) {
      return
    }
    try {
      await navigator.clipboard.writeText(enrollToken.token)
      setTokenCopied(true)
      window.setTimeout(() => setTokenCopied(false), 1500)
    } catch (_err) {
      setTokenCopied(false)
    }
  }

  const handleParse = () => {
    if (!selectedFile?.qeInput) {
      setParseError('QE input is missing. Re-import the .in file.')
      return
    }
    void runParse(selectedFile.qeInput, selectedFile.structureId)
  }

  const handleRun = async () => {
    if (!selectedFile?.qeInput) {
      setRunError('QE input is missing. Re-import the .in file.')
      return
    }
    if (!session) {
      setRunError('Sign in to run ZPE.')
      return
    }
    if (!activeTargetId) {
      setRunError('Select an active compute target before running.')
      return
    }
    if (!parseResult) {
      setRunError('Parsing is not complete yet.')
      return
    }
    if (mobileIndices.size === 0) {
      setRunError('Select at least one mobile atom to run ZPE.')
      return
    }
    setRunError(null)
    setIsSubmitting(true)
    setJobStatus(null)
    setJobResult(null)
    try {
      const response = await createZpeJob({
        content: selectedFile.qeInput,
        mobile_indices: Array.from(mobileIndices).sort((a, b) => a - b),
        use_environ: useEnviron,
        input_dir: inputDir.trim() ? inputDir.trim() : null,
        calc_mode: calcMode,
        structure_id: selectedFile.structureId ?? null,
      })
      setJobId(response.job_id)
    } catch (err) {
      setRunError(
        err instanceof Error ? err.message : 'ZPE job submission failed.',
      )
    } finally {
      setIsSubmitting(false)
    }
  }

  const handleDownload = async (kind: 'summary' | 'freqs') => {
    if (!session) {
      setRunError('Sign in to download results.')
      return
    }
    if (!jobId) {
      setRunError('No completed job available for download.')
      return
    }
    setRunError(null)
    try {
      const content = await downloadZpeFile(jobId, kind)
      const suffix = kind === 'summary' ? 'summary' : 'freqs'
      const extension = kind === 'summary' ? 'txt' : 'csv'
      downloadTextFile(
        content,
        `zpe-${suffix}-${jobId.slice(0, 8)}.${extension}`,
        kind === 'summary' ? 'text/plain' : 'text/csv',
      )
    } catch (err) {
      setRunError(err instanceof Error ? err.message : 'Download failed.')
    }
  }

  const handleToggleMobile = (index: number, fixed: boolean) => {
    if (!parseResult) {
      return
    }
    if (fixed) {
      return
    }
    setMobileIndices((prev) => {
      const next = new Set(prev)
      if (next.has(index)) {
        next.delete(index)
      } else {
        next.add(index)
      }
      return next
    })
  }

  const handleMolstarToggle = (index: number) => {
    if (!parseResult) {
      return
    }
    handleToggleMobile(index, fixedIndexSet.has(index))
  }

  const handleSelectAll = () => {
    if (!parseResult) {
      return
    }
    const next = new Set<number>()
    atomRows.forEach((row) => {
      if (!fixedIndexSet.has(row.index)) {
        next.add(row.index)
      }
    })
    setMobileIndices(next)
  }

  const handleSelectNone = () => {
    if (!parseResult) {
      return
    }
    setMobileIndices(new Set())
  }

  const handleInvertSelection = () => {
    if (!parseResult) {
      return
    }
    setMobileIndices((prev) => {
      const next = new Set<number>()
      atomRows.forEach((row) => {
        if (fixedIndexSet.has(row.index)) {
          return
        }
        if (!prev.has(row.index)) {
          next.add(row.index)
        }
      })
      return next
    })
  }

  const handleResetSelection = () => {
    if (!parseResult) {
      return
    }
    const fixedSet = new Set(parseResult.fixed_indices)
    const next = new Set<number>()
    parseResult.structure.atoms.forEach((_atom, index) => {
      if (!fixedSet.has(index)) {
        next.add(index)
      }
    })
    setMobileIndices(next)
  }

  const atomCount = atomRows.length
  const fixedCount = fixedIndexSet.size
  const mobileCount = mobileIndices.size
  const hasStructure = atomCount > 0
  const hasParseMeta = Boolean(parseResult)
  const selectionEnabled = hasStructure && hasParseMeta
  const selectionColorEnabled = selectionEnabled

  const canRun =
    Boolean(selectedFile?.qeInput) &&
    hasStructure &&
    hasParseMeta &&
    mobileIndices.size > 0 &&
    Boolean(session) &&
    Boolean(activeTargetId)

  return (
    <div className="grid min-h-0 flex-1 grid-cols-1 gap-4 overflow-hidden xl:grid-cols-[22rem_minmax(0,1fr)_22rem]">
      <div className="flex min-h-0 flex-col gap-3">
        <div className="flex min-h-0 flex-1 flex-col overflow-hidden rounded-lg border border-slate-200 bg-white shadow-sm">
          <div className="border-b border-slate-100 px-3 py-2">
            <div className="flex items-center justify-between">
              <div>
                <p className="text-xs uppercase tracking-[0.3em] text-slate-400">
                  3D Viewer
                </p>
                <p className="text-sm font-semibold text-slate-800">
                  Mol* (ASE Atoms)
                </p>
              </div>
              <Badge
                variant="outline"
                className="rounded-full border-slate-200 px-2 text-[10px] uppercase tracking-wide text-slate-500"
              >
                {viewerCifUrl ? 'cif' : 'missing'}
              </Badge>
            </div>
          </div>
          <div className="relative min-h-0 flex-1">
            {viewerCifUrl ? (
              <MolstarViewer
                cifUrl={viewerCifUrl}
                onError={setViewerError}
                onLoad={() => setViewerError(null)}
                selectedAtomIndices={hasParseMeta ? mobileIndexList : undefined}
                disabledAtomIndices={
                  hasParseMeta ? disabledAtomIndices : undefined
                }
                onAtomToggle={hasParseMeta ? handleMolstarToggle : undefined}
                className="h-full w-full rounded-none border-0"
              />
            ) : (
              <div className="flex h-full w-full flex-col items-center justify-center px-4 text-center text-muted-foreground">
                <Layers className="mb-2 h-12 w-12 text-slate-300" />
                <span className="font-medium text-slate-600">
                  No structure loaded
                </span>
                <span className="text-xs text-slate-400">
                  Import a QE input to fetch CIF.
                </span>
              </div>
            )}
            {viewerError ? (
              <div className="absolute inset-0 z-20 flex items-center justify-center bg-white/80 px-4 text-center">
                <div className="rounded-md border border-red-200 bg-red-50 px-3 py-2 text-sm text-red-700 shadow-sm">
                  <p className="font-semibold">Viewer failed to load</p>
                  <p className="mt-1 text-xs text-red-600">{viewerError}</p>
                </div>
              </div>
            ) : null}
          </div>
        </div>
      </div>

      <div className="flex min-h-0 flex-col gap-3 overflow-hidden">
        <div className="flex min-h-0 flex-1 flex-col overflow-hidden rounded-lg border border-slate-200 bg-white shadow-sm">
          <div className="border-b border-slate-100 px-3 py-2">
            <div className="flex flex-wrap items-center justify-between gap-2">
              <div>
                <p className="text-xs uppercase tracking-[0.3em] text-slate-400">
                  Atom Table
                </p>
                <p className="text-sm font-semibold text-slate-800">
                  Mobile Selection
                </p>
              </div>
              <div className="flex flex-wrap items-center gap-2">
                <Input
                  value={atomFilter}
                  onChange={(event) => setAtomFilter(event.target.value)}
                  placeholder="Filter by element or index"
                  className="h-8 w-44 text-xs"
                  disabled={!hasStructure}
                />
                <div className="flex flex-wrap gap-1">
                  <Button
                    variant="outline"
                    size="sm"
                    className="h-8 px-2 text-[11px]"
                    onClick={handleSelectAll}
                    disabled={!selectionEnabled}
                  >
                    All
                  </Button>
                  <Button
                    variant="outline"
                    size="sm"
                    className="h-8 px-2 text-[11px]"
                    onClick={handleSelectNone}
                    disabled={!selectionEnabled}
                  >
                    None
                  </Button>
                  <Button
                    variant="outline"
                    size="sm"
                    className="h-8 px-2 text-[11px]"
                    onClick={handleInvertSelection}
                    disabled={!selectionEnabled}
                  >
                    Invert
                  </Button>
                </div>
              </div>
            </div>
          </div>

          <div className="flex min-h-0 flex-1 flex-col">
            {hasStructure ? (
              <AtomTable
                rows={filteredRows}
                selectedIndices={
                  selectionColorEnabled ? mobileIndices : undefined
                }
                fixedIndices={selectionColorEnabled ? fixedIndexSet : undefined}
                onRowClick={
                  selectionEnabled
                    ? (index) =>
                        handleToggleMobile(index, fixedIndexSet.has(index))
                    : undefined
                }
                selectionEnabled={selectionEnabled}
                digits={3}
                emptyText="No matching atoms found."
                stickyHeader
                containerClassName="h-full max-h-full"
                showIndex
              />
            ) : (
              <div className="flex flex-1 items-center justify-center text-sm text-slate-500">
                Select a workspace file to load atoms.
              </div>
            )}
          </div>
        </div>
      </div>

      <div className="flex min-h-0 flex-col gap-3 overflow-y-auto pr-1">
        <div className="rounded-lg border border-slate-200 bg-white p-3 shadow-sm">
          <div className="flex items-center justify-between">
            <div>
              <p className="text-xs uppercase tracking-[0.3em] text-slate-400">
                Account
              </p>
              <p className="text-sm font-semibold text-slate-800">
                Sign in & Compute
              </p>
            </div>
            <Badge
              variant="outline"
              className="rounded-full border-slate-200 px-2 text-[10px] uppercase tracking-wide text-slate-500"
            >
              {session ? 'signed in' : 'guest'}
            </Badge>
          </div>
          <div className="mt-3 space-y-3 text-xs text-slate-600">
            {!session ? (
              <div className="space-y-2">
                <div className="flex gap-2">
                  <Button
                    type="button"
                    variant={authMode === 'login' ? 'default' : 'outline'}
                    size="sm"
                    className="h-7 px-2 text-[11px]"
                    onClick={() => {
                      setAuthMode('login')
                      setAuthError(null)
                    }}
                  >
                    Sign in
                  </Button>
                  <Button
                    type="button"
                    variant={authMode === 'register' ? 'default' : 'outline'}
                    size="sm"
                    className="h-7 px-2 text-[11px]"
                    onClick={() => {
                      setAuthMode('register')
                      setAuthError(null)
                    }}
                  >
                    Create
                  </Button>
                </div>
                <div className="grid gap-2">
                  <Input
                    value={authEmail}
                    onChange={(event) => setAuthEmail(event.target.value)}
                    placeholder="Email"
                    type="email"
                    className="h-8 text-xs"
                  />
                  <Input
                    value={authPassword}
                    onChange={(event) => setAuthPassword(event.target.value)}
                    placeholder="Password (min 8 chars)"
                    type="password"
                    className="h-8 text-xs"
                  />
                </div>
                <Button
                  type="button"
                  size="sm"
                  className="h-8 w-full text-xs font-semibold"
                  onClick={handleAuthSubmit}
                  disabled={authBusy}
                >
                  {authBusy
                    ? 'Working...'
                    : authMode === 'login'
                      ? 'Sign in'
                      : 'Create account'}
                </Button>
                {authError ? (
                  <div className="rounded-md border border-red-200 bg-red-50 px-3 py-2 text-[11px] text-red-600">
                    {authError}
                  </div>
                ) : null}
              </div>
            ) : (
              <div className="space-y-2">
                <div className="rounded-md border border-slate-100 bg-slate-50 px-3 py-2 text-[11px] text-slate-600">
                  <p className="text-[10px] uppercase tracking-[0.3em] text-slate-400">
                    Signed in
                  </p>
                  <p className="mt-1 font-medium text-slate-700">
                    {session.user.email}
                  </p>
                  <p className="mt-1 text-[10px] text-slate-400">
                    expires {session.expires_at}
                  </p>
                </div>
                <div className="flex flex-wrap gap-2">
                  <Button
                    type="button"
                    variant="outline"
                    size="sm"
                    className="h-7 px-2 text-[11px]"
                    onClick={refreshTargets}
                    disabled={targetsBusy}
                  >
                    Refresh targets
                  </Button>
                  <Button
                    type="button"
                    variant="outline"
                    size="sm"
                    className="h-7 px-2 text-[11px]"
                    onClick={handleLogout}
                    disabled={authBusy}
                  >
                    Sign out
                  </Button>
                </div>
              </div>
            )}

            <div className="rounded-md border border-slate-100 bg-slate-50 px-3 py-2">
              <p className="text-[10px] uppercase tracking-[0.3em] text-slate-400">
                Compute queue
              </p>
              {!session ? (
                <p className="mt-2 text-[11px] text-slate-400">
                  Sign in to register and select your compute targets.
                </p>
              ) : (
                <div className="mt-2 space-y-2">
                  <div>
                    <label className="text-[11px] font-medium text-slate-500">
                      Active target
                    </label>
                    <Select
                      value={activeTargetId ?? ''}
                      onValueChange={handleSelectTarget}
                      disabled={queueTargets.length === 0 || targetsBusy}
                    >
                      <SelectTrigger className="mt-1 h-8 text-xs">
                        <SelectValue
                          placeholder={
                            queueTargets.length > 0
                              ? 'Select target'
                              : 'No targets yet'
                          }
                        />
                      </SelectTrigger>
                      <SelectContent>
                        {queueTargets.map((target) => (
                          <SelectItem
                            key={target.target_id}
                            value={target.target_id}
                          >
                            {target.name ?? target.queue_name} ·{' '}
                            {target.queue_name}
                          </SelectItem>
                        ))}
                      </SelectContent>
                    </Select>
                  </div>

                  {activeTarget ? (
                    <div className="text-[11px] text-slate-500">
                      Queue:{' '}
                      <span className="font-mono">
                        {activeTarget.queue_name}
                      </span>
                    </div>
                  ) : (
                    <p className="text-[11px] text-slate-400">
                      Register a compute worker to activate a queue.
                    </p>
                  )}

                  <div className="flex flex-wrap items-center gap-2">
                    <Button
                      type="button"
                      size="sm"
                      className="h-7 px-2 text-[11px]"
                      onClick={handleGenerateToken}
                      disabled={enrollBusy}
                    >
                      {enrollBusy ? 'Generating...' : 'Generate enroll token'}
                    </Button>
                    {enrollToken ? (
                      <Button
                        type="button"
                        variant="outline"
                        size="sm"
                        className="h-7 px-2 text-[11px]"
                        onClick={handleCopyToken}
                      >
                        <Copy className="h-3 w-3" />
                        {tokenCopied ? 'Copied' : 'Copy'}
                      </Button>
                    ) : null}
                  </div>

                  {enrollToken ? (
                    <div className="rounded-md border border-slate-200 bg-white px-2 py-2 text-[10px] text-slate-500">
                      <p className="text-[10px] uppercase tracking-[0.3em] text-slate-400">
                        Enroll token
                      </p>
                      <p className="mt-1 break-all font-mono text-[10px] text-slate-700">
                        {enrollToken.token}
                      </p>
                      <p className="mt-1 text-[10px] text-slate-400">
                        expires {enrollToken.expires_at}
                      </p>
                    </div>
                  ) : null}

                  {enrollError ? (
                    <div className="rounded-md border border-red-200 bg-red-50 px-3 py-2 text-[11px] text-red-600">
                      {enrollError}
                    </div>
                  ) : null}

                  {targetsError ? (
                    <div className="rounded-md border border-amber-200 bg-amber-50 px-3 py-2 text-[11px] text-amber-700">
                      {targetsError}
                    </div>
                  ) : null}
                </div>
              )}
            </div>
          </div>
        </div>

        <div className="rounded-lg border border-slate-200 bg-white p-3 shadow-sm">
          <div className="flex items-center justify-between">
            <div>
              <p className="text-xs uppercase tracking-[0.3em] text-slate-400">
                ZPE Input
              </p>
              <p className="text-sm font-semibold text-slate-800">
                QE .in Source
              </p>
            </div>
            <Badge
              variant="outline"
              className="rounded-full border-slate-200 px-2 text-[10px] uppercase tracking-wide text-slate-500"
            >
              {selectedFile?.qeInput ? 'ready' : 'missing'}
            </Badge>
          </div>
          <div className="mt-3 space-y-3">
            <div>
              <label className="text-xs font-medium text-slate-500">
                Workspace file
              </label>
              <Select
                value={selectedFile?.id}
                onValueChange={(value) => setSelectedFileId(value)}
                disabled={availableFiles.length === 0}
              >
                <SelectTrigger className="mt-1">
                  <SelectValue
                    placeholder={
                      availableFiles.length > 0
                        ? 'Select .in file'
                        : 'No .in file'
                    }
                  />
                </SelectTrigger>
                <SelectContent>
                  {availableFiles.map((file) => (
                    <SelectItem key={file.id} value={file.id}>
                      {file.name}
                    </SelectItem>
                  ))}
                </SelectContent>
              </Select>
            </div>

            <div className="flex flex-wrap items-center gap-2">
              <Button
                size="sm"
                className="h-8 bg-rose-600 px-3 text-xs font-semibold text-white hover:bg-rose-700"
                onClick={handleParse}
                disabled={isParsing || !selectedFile?.qeInput}
              >
                {isParsing ? 'Parsing...' : 'Refresh Parse'}
              </Button>
              <Button
                variant="outline"
                size="sm"
                className="h-8 px-3 text-xs"
                onClick={handleResetSelection}
                disabled={!selectionEnabled}
              >
                <RefreshCw className="h-3 w-3" />
                Reset Mobile
              </Button>
            </div>

            {parseError ? (
              <div className="rounded-md border border-red-200 bg-red-50 px-3 py-2 text-xs text-red-600">
                {parseError}
              </div>
            ) : null}
          </div>
        </div>

        <div className="rounded-lg border border-slate-200 bg-white p-3 shadow-sm">
          <div className="flex items-center justify-between">
            <div>
              <p className="text-xs uppercase tracking-[0.3em] text-slate-400">
                Parse Summary
              </p>
              <p className="text-sm font-semibold text-slate-800">
                Fixed / Mobile Atoms
              </p>
            </div>
            <div className="flex items-center gap-1 text-xs text-slate-500">
              <span>{atomCount} atoms</span>
            </div>
          </div>
          <div className="mt-3 space-y-3 text-xs text-slate-600">
            <div className="flex items-center justify-between rounded-md border border-slate-100 bg-slate-50 px-3 py-2">
              <span>Fixed</span>
              <span className="font-semibold">{fixedCount}</span>
            </div>
            <div className="flex items-center justify-between rounded-md border border-slate-100 bg-slate-50 px-3 py-2">
              <span>Mobile</span>
              <span className="font-semibold">{mobileCount}</span>
            </div>
            <div className="flex items-center justify-between rounded-md border border-slate-100 bg-slate-50 px-3 py-2">
              <span>K-points</span>
              <span className="font-mono text-[11px]">
                {hasParseMeta
                  ? parseResult?.kpoints
                    ? parseResult.kpoints.join(' × ')
                    : 'n/a'
                  : isParsing
                    ? 'parsing…'
                    : '—'}
              </span>
            </div>
            <div>
              <p className="text-[10px] uppercase tracking-[0.3em] text-slate-400">
                Atomic species
              </p>
              <div className="mt-2 flex flex-wrap gap-1">
                {!hasParseMeta ? (
                  <span className="text-[11px] text-slate-400">
                    {isParsing ? 'Parsing…' : '—'}
                  </span>
                ) : speciesEntries.length > 0 ? (
                  speciesEntries.map(([symbol, pseudo]) => (
                    <span
                      key={symbol}
                      className="rounded-full border border-slate-200 bg-white px-2 py-0.5 text-[10px] text-slate-500"
                      title={pseudo}
                    >
                      {symbol}: {pseudo}
                    </span>
                  ))
                ) : (
                  <span className="text-[11px] text-slate-400">
                    ATOMIC_SPECIES not found
                  </span>
                )}
              </div>
            </div>
          </div>
        </div>

        <div className="rounded-lg border border-slate-200 bg-white p-3 text-xs text-slate-500 shadow-sm">
          <div className="flex items-center gap-2">
            <Activity className="h-4 w-4 text-rose-500" />
            <span className="font-semibold text-slate-700">
              ZPE selection notes
            </span>
          </div>
          <ul className="mt-2 space-y-1 text-[11px] text-slate-500">
            <li>Fixed atoms are locked from mobile selection.</li>
            <li>Selection is 1-based in UI, 0-based in API.</li>
            <li>Parse to refresh mobile defaults from flags.</li>
          </ul>
        </div>

        <div className="rounded-lg border border-slate-200 bg-white p-3 shadow-sm">
          <div className="flex items-center justify-between">
            <div>
              <p className="text-xs uppercase tracking-[0.3em] text-slate-400">
                Run
              </p>
              <p className="text-sm font-semibold text-slate-800">
                ZPE Job Settings
              </p>
            </div>
          </div>
          <div className="mt-3 space-y-3 text-xs text-slate-600">
            <div className="flex items-center justify-between rounded-md border border-slate-100 bg-slate-50 px-3 py-2">
              <span>Use environ.in</span>
              <Switch checked={useEnviron} onCheckedChange={setUseEnviron} />
            </div>
            <div>
              <label className="text-xs font-medium text-slate-500">
                Calc mode
              </label>
              <Select
                value={calcMode}
                onValueChange={(value) =>
                  setCalcMode(value as 'new' | 'continue')
                }
              >
                <SelectTrigger className="mt-1 h-8 text-xs">
                  <SelectValue />
                </SelectTrigger>
                <SelectContent>
                  <SelectItem value="continue">continue</SelectItem>
                  <SelectItem value="new">new</SelectItem>
                </SelectContent>
              </Select>
            </div>

            <CollapsibleSection title="Advanced" defaultOpen={false}>
              <div className="space-y-2">
                <label className="text-xs font-medium text-slate-500">
                  Input directory
                </label>
                <Input
                  value={inputDir}
                  onChange={(event) => setInputDir(event.target.value)}
                  placeholder="Optional input_dir"
                  className="h-8 text-xs"
                />
              </div>
            </CollapsibleSection>

            <Button
              className="mt-2 h-10 w-full bg-rose-600 text-sm font-semibold text-white hover:bg-rose-700"
              onClick={handleRun}
              disabled={!canRun || isSubmitting}
            >
              <Play className="h-4 w-4" />
              {isSubmitting ? 'Submitting...' : 'Run ZPE'}
            </Button>
            {!session ? (
              <p className="text-[11px] text-slate-400">
                Sign in to enable remote compute.
              </p>
            ) : !activeTargetId ? (
              <p className="text-[11px] text-slate-400">
                Select a compute target before running.
              </p>
            ) : null}
            {!selectedFile?.qeInput ? (
              <p className="text-[11px] text-slate-400">
                Import a QE .in file to enable ZPE.
              </p>
            ) : null}
          </div>
        </div>

        <div className="rounded-lg border border-slate-200 bg-white p-3 shadow-sm">
          <div className="flex items-center justify-between">
            <div>
              <p className="text-xs uppercase tracking-[0.3em] text-slate-400">
                Status
              </p>
              <p className="text-sm font-semibold text-slate-800">
                Job Monitor
              </p>
            </div>
            <span
              className={cn(
                'inline-flex items-center gap-1 rounded-full border px-2 py-0.5 text-[10px] font-semibold uppercase tracking-wide',
                statusTone(jobStatus?.status),
              )}
            >
              {jobStatus?.status === 'finished' ? (
                <CheckCircle2 className="h-3 w-3" />
              ) : jobStatus?.status === 'failed' ? (
                <AlertTriangle className="h-3 w-3" />
              ) : null}
              {statusLabel(jobStatus?.status)}
            </span>
          </div>
          <div className="mt-3 space-y-2 text-xs text-slate-600">
            <div className="rounded-md border border-slate-100 bg-slate-50 px-3 py-2">
              <p className="text-[10px] uppercase tracking-[0.3em] text-slate-400">
                Job ID
              </p>
              <p className="mt-1 break-all font-mono text-[11px]">
                {jobId ?? '—'}
              </p>
            </div>
            <div className="rounded-md border border-slate-100 bg-slate-50 px-3 py-2">
              <p className="text-[10px] uppercase tracking-[0.3em] text-slate-400">
                Updated
              </p>
              <p className="mt-1 text-[11px]">{jobStatus?.updated_at ?? '—'}</p>
            </div>
            {jobStatus?.detail ? (
              <div className="rounded-md border border-amber-200 bg-amber-50 px-3 py-2 text-[11px] text-amber-700">
                {jobStatus.detail}
              </div>
            ) : null}
            {runError ? (
              <div className="rounded-md border border-red-200 bg-red-50 px-3 py-2 text-[11px] text-red-600">
                {runError}
              </div>
            ) : null}
          </div>
        </div>

        <div className="rounded-lg border border-slate-200 bg-white p-3 shadow-sm">
          <div className="flex items-center justify-between">
            <div>
              <p className="text-xs uppercase tracking-[0.3em] text-slate-400">
                Result
              </p>
              <p className="text-sm font-semibold text-slate-800">
                ZPE Summary
              </p>
            </div>
          </div>
          <div className="mt-3 space-y-3">
            {jobResult ? (
              <div className="grid gap-2 text-xs text-slate-600">
                <div className="grid grid-cols-2 gap-2">
                  <div className="rounded-md border border-slate-100 bg-slate-50 px-3 py-2">
                    <p className="text-[10px] uppercase tracking-[0.3em] text-slate-400">
                      ZPE (eV)
                    </p>
                    <p className="mt-1 text-sm font-semibold text-slate-700">
                      {formatNumber(jobResult.zpe_ev, 6)}
                    </p>
                  </div>
                  <div className="rounded-md border border-slate-100 bg-slate-50 px-3 py-2">
                    <p className="text-[10px] uppercase tracking-[0.3em] text-slate-400">
                      S_vib
                    </p>
                    <p className="mt-1 text-sm font-semibold text-slate-700">
                      {formatNumber(jobResult.s_vib_jmol_k, 3)}
                    </p>
                  </div>
                </div>

                <div className="grid grid-cols-2 gap-2">
                  <div className="rounded-md border border-slate-100 bg-slate-50 px-3 py-2">
                    <p className="text-[10px] uppercase tracking-[0.3em] text-slate-400">
                      Elapsed
                    </p>
                    <p className="mt-1 text-sm font-semibold text-slate-700">
                      {formatNumber(jobResult.elapsed_seconds, 2)}s
                    </p>
                  </div>
                  <div className="rounded-md border border-slate-100 bg-slate-50 px-3 py-2">
                    <p className="text-[10px] uppercase tracking-[0.3em] text-slate-400">
                      Modes
                    </p>
                    <p className="mt-1 text-sm font-semibold text-slate-700">
                      {jobResult.freqs_cm.length}
                    </p>
                  </div>
                </div>

                <div className="rounded-md border border-slate-100 bg-slate-50 px-3 py-2 text-[11px] text-slate-500">
                  <div className="flex items-center justify-between">
                    <span>freqs min / max</span>
                    <span className="font-mono">
                      {freqSummary
                        ? `${formatNumber(freqSummary.min, 2)} / ${formatNumber(freqSummary.max, 2)}`
                        : '—'}
                    </span>
                  </div>
                  <div className="mt-1 flex items-center justify-between">
                    <span>imaginary modes</span>
                    <span className="font-mono">
                      {freqSummary ? freqSummary.imagCount : '—'}
                    </span>
                  </div>
                </div>
              </div>
            ) : (
              <div className="rounded-md border border-dashed border-slate-200 bg-slate-50 px-3 py-4 text-center text-xs text-slate-500">
                Results appear after the job finishes.
              </div>
            )}

            <div className="flex flex-wrap gap-2">
              <Button
                variant="outline"
                size="sm"
                className="h-8 text-xs"
                onClick={() => handleDownload('summary')}
                disabled={!jobResult}
              >
                <Download className="h-3 w-3" />
                summary.txt
              </Button>
              <Button
                variant="outline"
                size="sm"
                className="h-8 text-xs"
                onClick={() => handleDownload('freqs')}
                disabled={!jobResult}
              >
                <Download className="h-3 w-3" />
                freqs.csv
              </Button>
            </div>
          </div>
        </div>
      </div>
    </div>
  )
}

function TransferToolPanel({
  structures = [],
}: {
  structures?: Array<WorkspaceFile>
}) {
  const availableStructures = useMemo(
    () => structures.filter((file) => file.structure || file.structureId),
    [structures],
  )
  const [sourceId, setSourceId] = useState<string | null>(null)
  const [targetId, setTargetId] = useState<string | null>(null)
  const [transferSummary, setTransferSummary] =
    useState<TransferSummary | null>(null)
  const [exportContent, setExportContent] = useState<string | null>(null)
  const [exportFilename, setExportFilename] = useState('')
  const [transferError, setTransferError] = useState<string | null>(null)
  const [viewerError, setViewerError] = useState<string | null>(null)
  const [isApplying, setIsApplying] = useState(false)
  const structureCacheRef = useRef<Partial<Record<string, Structure>>>({})
  const applyTokenRef = useRef(0)

  const fileById = useMemo(
    () => new Map(availableStructures.map((file) => [file.id, file])),
    [availableStructures],
  )
  const sourceFile = sourceId ? (fileById.get(sourceId) ?? null) : null
  const targetFile = targetId ? (fileById.get(targetId) ?? null) : null

  useEffect(() => {
    if (availableStructures.length === 0) {
      setSourceId(null)
      setTargetId(null)
      return
    }
    let nextSourceId = sourceId
    if (!nextSourceId || !fileById.has(nextSourceId)) {
      nextSourceId = availableStructures[0].id
      if (nextSourceId !== sourceId) {
        setSourceId(nextSourceId)
      }
    }
    if (!targetId || !fileById.has(targetId) || targetId === nextSourceId) {
      const fallback = availableStructures.find(
        (file) => file.id !== nextSourceId,
      )
      const nextTargetId = fallback?.id ?? null
      if (nextTargetId !== targetId) {
        setTargetId(nextTargetId)
      }
    }
  }, [availableStructures, fileById, sourceId, targetId])

  useEffect(() => {
    setTransferSummary(null)
    setExportContent(null)
    setExportFilename('')
    setTransferError(null)
  }, [sourceId, targetId])

  useEffect(() => {
    setViewerError(null)
  }, [transferSummary?.pdbText])

  const resolveStructure = useCallback(async (file: WorkspaceFile) => {
    if (file.structure) {
      return file.structure
    }
    const cached = structureCacheRef.current[file.id]
    if (cached) {
      return cached
    }
    if (!file.structureId) {
      return null
    }
    const nextStructure = await getStructure(file.structureId)
    structureCacheRef.current[file.id] = nextStructure
    return nextStructure
  }, [])

  const handleApply = useCallback(async () => {
    if (!sourceId || !targetId) {
      setTransferError('Select both source and target structures.')
      return
    }
    if (sourceId === targetId) {
      setTransferError('Source and target must be different.')
      return
    }
    if (!sourceFile || !targetFile) {
      setTransferError('Selected structures are unavailable.')
      return
    }
    setTransferError(null)
    setExportContent(null)
    setExportFilename('')
    const token = applyTokenRef.current + 1
    applyTokenRef.current = token
    setIsApplying(true)
    try {
      const [sourceStructure, targetStructure] = await Promise.all([
        resolveStructure(sourceFile),
        resolveStructure(targetFile),
      ])
      if (applyTokenRef.current !== token) {
        return
      }
      if (!sourceStructure || !targetStructure) {
        setTransferError('Missing structure data for transfer.')
        return
      }
      const sourceAtoms = sourceStructure.atoms.length
      const targetAtoms = targetStructure.atoms.length
      if (sourceAtoms === 0 || targetAtoms === 0) {
        setTransferError('Source and target must contain atoms.')
        return
      }
      const transferredAtoms = Math.min(sourceAtoms, targetAtoms)
      const nextAtoms = targetStructure.atoms.map((atom, index) => {
        if (index >= transferredAtoms) {
          return atom
        }
        const sourceAtom = sourceStructure.atoms[index]
        return {
          ...atom,
          x: sourceAtom.x,
          y: sourceAtom.y,
          z: sourceAtom.z,
        }
      })
      const nextStructure: Structure = {
        ...targetStructure,
        atoms: nextAtoms,
      }
      const pdbText = atomsToPdb(nextAtoms)
      setTransferSummary({
        structure: nextStructure,
        pdbText,
        sourceAtoms,
        targetAtoms,
        transferredAtoms,
      })
      const content = await exportQeInput(nextStructure)
      if (applyTokenRef.current !== token) {
        return
      }
      setExportContent(content)
      setExportFilename(createTransferFilename(targetFile.name))
    } catch (err) {
      if (applyTokenRef.current !== token) {
        return
      }
      setTransferError(
        err instanceof Error ? err.message : 'Transfer failed to run.',
      )
    } finally {
      if (applyTokenRef.current === token) {
        setIsApplying(false)
      }
    }
  }, [resolveStructure, sourceFile, sourceId, targetFile, targetId])

  const handleDownload = () => {
    if (!exportContent) {
      setTransferError('No transferred .in file is ready to download.')
      return
    }
    const filename = exportFilename || 'Transferred.in'
    downloadTextFile(exportContent, filename, 'text/plain')
  }

  const previewReady = Boolean(transferSummary?.pdbText)
  const canApply =
    Boolean(sourceId && targetId && sourceId !== targetId) &&
    !isApplying &&
    availableStructures.length > 1
  const canDownload = Boolean(exportContent)
  const summarySource =
    transferSummary?.sourceAtoms ?? sourceFile?.structure?.atoms.length ?? null
  const summaryTarget =
    transferSummary?.targetAtoms ?? targetFile?.structure?.atoms.length ?? null
  const summaryTransferred = transferSummary?.transferredAtoms ?? null

  return (
    <div className="flex flex-1 gap-4 overflow-hidden">
      <div className="flex w-[26rem] flex-shrink-0 flex-col gap-4 overflow-y-auto pr-1">
        <div className="flex min-h-[260px] flex-1 flex-col overflow-hidden rounded-lg border border-slate-200 bg-white shadow-sm">
          <div className="border-b border-slate-100 px-3 py-2">
            <div className="flex items-center justify-between">
              <div>
                <p className="text-xs uppercase tracking-[0.3em] text-slate-400">
                  Action Preview
                </p>
                <p className="text-sm font-semibold text-slate-800">
                  Transferred Structure
                </p>
              </div>
              <Badge
                variant="outline"
                className="rounded-full border-slate-200 px-2 text-[10px] uppercase tracking-wide text-slate-500"
              >
                {previewReady ? 'ready' : 'idle'}
              </Badge>
            </div>
          </div>
          <div className="relative min-h-0 flex-1">
            {previewReady && transferSummary ? (
              <MolstarViewer
                pdbText={transferSummary.pdbText}
                onError={setViewerError}
                onLoad={() => setViewerError(null)}
                className="h-full w-full rounded-none border-0"
              />
            ) : (
              <div className="flex h-full w-full flex-col items-center justify-center px-4 text-center text-muted-foreground">
                <Layers className="mb-2 h-12 w-12 text-slate-300" />
                <span className="font-medium text-slate-600">
                  No transfer preview yet
                </span>
                <span className="text-xs text-slate-400">
                  Select source/target and apply transfer.
                </span>
              </div>
            )}
            {viewerError ? (
              <div className="absolute inset-0 z-20 flex items-center justify-center bg-white/80 px-4 text-center">
                <div className="rounded-md border border-red-200 bg-red-50 px-3 py-2 text-sm text-red-700 shadow-sm">
                  <p className="font-semibold">Viewer failed to load</p>
                  <p className="mt-1 text-xs text-red-600">{viewerError}</p>
                </div>
              </div>
            ) : null}
          </div>
        </div>

        <div className="rounded-lg border border-slate-200 bg-white p-3 shadow-sm">
          <div className="flex items-center justify-between">
            <div>
              <p className="text-xs uppercase tracking-[0.3em] text-slate-400">
                Transfer Summary
              </p>
              <p className="text-sm font-semibold text-slate-800">
                Result Overview
              </p>
            </div>
            <Badge
              variant="outline"
              className="rounded-full border-slate-200 px-2 text-[10px] uppercase tracking-wide text-slate-500"
            >
              {transferSummary ? 'ready' : 'pending'}
            </Badge>
          </div>
          <div className="mt-3 space-y-2 text-xs text-slate-600">
            <div className="flex items-center justify-between rounded-md border border-slate-100 bg-slate-50 px-3 py-2">
              <span>Source atoms</span>
              <span className="font-mono">{summarySource ?? '—'}</span>
            </div>
            <div className="flex items-center justify-between rounded-md border border-slate-100 bg-slate-50 px-3 py-2">
              <span>Target atoms</span>
              <span className="font-mono">{summaryTarget ?? '—'}</span>
            </div>
            <div className="flex items-center justify-between rounded-md border border-slate-100 bg-slate-50 px-3 py-2">
              <span>Transferred</span>
              <span className="font-mono">{summaryTransferred ?? '—'}</span>
            </div>
          </div>
          <div className="mt-3 flex flex-wrap gap-2">
            <Button
              variant="outline"
              size="sm"
              className="h-8 text-xs"
              onClick={handleDownload}
              disabled={!canDownload}
            >
              <Download className="h-3 w-3" />
              Download .in
            </Button>
          </div>
        </div>
      </div>

      <div className="flex min-w-0 flex-1 flex-col overflow-y-auto">
        <div className="flex h-full flex-col overflow-hidden rounded-md border border-blue-200 bg-white shadow-sm">
          <div className="flex shrink-0 items-center gap-2 border-b border-blue-100 bg-blue-50/50 px-3 py-2">
            <MousePointerClick className="h-4 w-4 text-blue-600" />
            <span className="text-sm font-medium text-blue-900">Actions</span>
          </div>
          <div className="flex-1 overflow-y-auto p-3">
            <div className="space-y-4">
              <div>
                <label className="text-xs font-medium text-slate-500">
                  Source structure
                </label>
                <Select
                  value={sourceId ?? undefined}
                  onValueChange={(value) => setSourceId(value)}
                  disabled={availableStructures.length === 0}
                >
                  <SelectTrigger className="mt-1">
                    <SelectValue
                      placeholder={
                        availableStructures.length > 0
                          ? 'Select source'
                          : 'No structures'
                      }
                    />
                  </SelectTrigger>
                  <SelectContent>
                    {availableStructures.map((file) => (
                      <SelectItem
                        key={file.id}
                        value={file.id}
                        disabled={file.id === targetId}
                      >
                        <div className="flex w-full items-center justify-between gap-2">
                          <span className="truncate">{file.name}</span>
                          <span className="rounded bg-slate-100 px-1.5 py-0.5 text-[10px] uppercase text-slate-500">
                            {file.kind}
                          </span>
                        </div>
                      </SelectItem>
                    ))}
                  </SelectContent>
                </Select>
              </div>

              <div>
                <label className="text-xs font-medium text-slate-500">
                  Target structure
                </label>
                <Select
                  value={targetId ?? undefined}
                  onValueChange={(value) => setTargetId(value)}
                  disabled={availableStructures.length < 2}
                >
                  <SelectTrigger className="mt-1">
                    <SelectValue
                      placeholder={
                        availableStructures.length > 1
                          ? 'Select target'
                          : 'Need at least two structures'
                      }
                    />
                  </SelectTrigger>
                  <SelectContent>
                    {availableStructures.map((file) => (
                      <SelectItem
                        key={file.id}
                        value={file.id}
                        disabled={file.id === sourceId}
                      >
                        <div className="flex w-full items-center justify-between gap-2">
                          <span className="truncate">{file.name}</span>
                          <span className="rounded bg-slate-100 px-1.5 py-0.5 text-[10px] uppercase text-slate-500">
                            {file.kind}
                          </span>
                        </div>
                      </SelectItem>
                    ))}
                  </SelectContent>
                </Select>
              </div>

              <div className="rounded-md border border-slate-100 bg-slate-50 px-3 py-2 text-xs text-slate-600">
                <div className="flex items-center justify-between">
                  <span>Source atoms</span>
                  <span className="font-mono">{summarySource ?? '—'}</span>
                </div>
                <div className="mt-1 flex items-center justify-between">
                  <span>Target atoms</span>
                  <span className="font-mono">{summaryTarget ?? '—'}</span>
                </div>
              </div>

              {transferError ? (
                <div className="rounded-md border border-red-200 bg-red-50 px-3 py-2 text-xs text-red-600">
                  <div className="flex items-center gap-2">
                    <AlertTriangle className="h-4 w-4" />
                    <span>{transferError}</span>
                  </div>
                </div>
              ) : null}

              {transferSummary && exportContent && !transferError ? (
                <div className="rounded-md border border-emerald-200 bg-emerald-50 px-3 py-2 text-xs text-emerald-700">
                  <div className="flex items-center gap-2">
                    <CheckCircle2 className="h-4 w-4" />
                    <span>
                      Transfer complete. Export ready for{' '}
                      {exportFilename || 'download'}.
                    </span>
                  </div>
                </div>
              ) : null}

              <div className="flex flex-wrap gap-2 pt-1">
                <Button
                  className="h-9 bg-emerald-600 px-4 text-xs font-semibold text-white hover:bg-emerald-700"
                  onClick={handleApply}
                  disabled={!canApply}
                >
                  {isApplying ? 'Applying…' : 'Apply Transfer'}
                </Button>
                <Button
                  variant="outline"
                  size="sm"
                  className="h-9 px-3 text-xs"
                  onClick={handleDownload}
                  disabled={!canDownload}
                >
                  <Download className="h-3 w-3" />
                  Download .in
                </Button>
              </div>
            </div>
          </div>
        </div>
      </div>
    </div>
  )
}

export function ToolPanel({
  mode,
  files,
  onClose,
  variant = 'stack',
  showHeader = true,
  showClose = true,
  className,
  structures,
  onSupercellCreated,
}: ToolPanelProps) {
  return (
    <div
      className={cn(
        'flex h-full flex-col gap-4 bg-slate-50/80 p-4',
        variant === 'stack'
          ? 'w-[60rem] flex-shrink-0 border-r border-border border-l-4 border-l-blue-500/20 animate-in fade-in slide-in-from-right-4 duration-300'
          : 'w-full',
        className,
      )}
    >
      {showHeader ? (
        <div className="flex shrink-0 items-center justify-between">
          <h2 className="text-lg font-semibold tracking-tight text-slate-800">
            {toolTitles[mode]}
          </h2>
          {showClose && onClose ? (
            <button
              type="button"
              onClick={onClose}
              className="rounded p-1 text-slate-400 transition-colors hover:bg-slate-200 hover:text-slate-600"
            >
              <X className="h-4 w-4" />
            </button>
          ) : null}
        </div>
      ) : null}

      {(() => {
        switch (mode) {
          case 'vibration':
            return <ZpeToolPanel files={files} />
          case 'supercell':
            return (
              <SupercellTool
                structures={structures ?? []}
                onSupercellCreated={onSupercellCreated}
              />
            )
          case 'transfer':
            return <TransferToolPanel structures={structures ?? files ?? []} />
          default:
            return null
        }
      })()}
    </div>
  )
}
