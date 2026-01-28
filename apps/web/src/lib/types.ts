export type Atom = {
  symbol: string
  x: number
  y: number
  z: number
}

export type Structure = {
  atoms: Array<Atom>
  lattice?: Lattice | null
}

export type Vector3 = {
  x: number
  y: number
  z: number
}

export type Lattice = {
  a: Vector3
  b: Vector3
  c: Vector3
}

export type LatticeParams = {
  a: number
  b: number
  c: number
  alpha: number
  beta: number
  gamma: number
}

export type SupercellMeta = {
  na: number
  nb: number
  layers: number
  overlapCount?: number
}
export type ZPEParseResponse = {
  structure: Structure
  fixed_indices: Array<number>
  atomic_species: Record<string, string>
  kpoints?: [number, number, number] | null
}

export type ZPEJobRequest = {
  content: string
  mobile_indices: Array<number>
  use_environ?: boolean
  input_dir?: string | null
  calc_mode?: 'new' | 'continue'
  structure_id?: string | null
}

export type ZPEJobResponse = {
  job_id: string
}

export type ZPEJobStatus = {
  status: string
  detail?: string | null
  updated_at?: string | null
}

export type AuthUser = {
  user_id: string
  email: string
  created_at: string
}

export type AuthSession = {
  token: string
  expires_at: string
  user: AuthUser
}

export type AuthMe = {
  user: AuthUser
  expires_at: string
}

export type ZPEQueueTarget = {
  target_id: string
  queue_name: string
  server_id: string
  registered_at: string
  name?: string | null
}

export type ZPEQueueTargetList = {
  targets: Array<ZPEQueueTarget>
  active_target_id?: string | null
}

export type ZPEResult = {
  freqs_cm: Array<number>
  zpe_ev: number
  s_vib_jmol_k: number
  mobile_indices: Array<number>
  fixed_indices: Array<number>
  kpts: [number, number, number]
  delta: number
  low_cut_cm: number
  temperature: number
  use_environ: boolean
  qe_input: string
  pseudo_dir: string
  calc_start_time: string
  calc_end_time: string
  elapsed_seconds: number
  cache_checked: number
  cache_deleted: number
  ecutwfc?: number | null
  ecutrho?: number | null
}
export type {
  SupercellGridAxis,
  SupercellGrid,
  SupercellBuildOptions,
  SupercellBuildOutput,
  SupercellBuildRequest,
  SupercellBuildMeta,
  SupercellBuildResponse,
} from '@chem-model/shared/types'
