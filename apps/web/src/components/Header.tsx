import { Link, useRouterState } from '@tanstack/react-router'

import { Activity, Database, Share2, Upload } from 'lucide-react'

import { Button } from '@/components/ui/button'

export default function Header() {
  const pathname = useRouterState({
    select: (state) => state.location.pathname,
  })

  if (pathname.startsWith('/editor-v2')) {
    return null
  }

  return (
    <header className="sticky top-0 z-40 border-b border-white/10 bg-slate-950/80 text-white backdrop-blur">
      <div className="mx-auto flex h-16 max-w-[1400px] items-center justify-between px-6">
        <div className="flex items-center gap-6">
          <Link to="/editor" className="flex items-center gap-3">
            <span className="flex h-9 w-9 items-center justify-center rounded-xl bg-gradient-to-br from-amber-200 via-orange-300 to-rose-400 text-slate-900">
              <Activity className="h-5 w-5" />
            </span>
            <div className="leading-tight">
              <p className="text-sm uppercase tracking-[0.3em] text-amber-200/70">
                Chem Lab
              </p>
              <p className="text-lg font-semibold text-white">
                Structure Studio
              </p>
            </div>
          </Link>
          <nav className="hidden items-center gap-4 text-sm text-white/70 md:flex">
            <Link
              to="/editor"
              className="rounded-full px-4 py-1.5 transition hover:text-white"
              activeProps={{
                className: 'rounded-full bg-white/10 px-4 py-1.5 text-white',
              }}
            >
              Editor
            </Link>
            <Link
              to="/editor-v2"
              className="rounded-full px-4 py-1.5 transition hover:text-white"
              activeProps={{
                className: 'rounded-full bg-white/10 px-4 py-1.5 text-white',
              }}
            >
              Editor V2
            </Link>
            <Link
              to="/supercell"
              className="rounded-full px-4 py-1.5 transition hover:text-white"
              activeProps={{
                className: 'rounded-full bg-white/10 px-4 py-1.5 text-white',
              }}
            >
              Supercell
            </Link>
            <Link
              to="/transplant"
              className="rounded-full px-4 py-1.5 transition hover:text-white"
              activeProps={{
                className: 'rounded-full bg-white/10 px-4 py-1.5 text-white',
              }}
            >
              Transplant
            </Link>
          </nav>
        </div>
        <div className="flex items-center gap-2 text-sm">
          <Button
            variant="outline"
            className="rounded-full border-white/15 bg-white/5 text-white/80 hover:bg-white/10"
            onClick={() =>
              window.dispatchEvent(new CustomEvent('chem-model-import'))
            }
          >
            <Upload className="h-4 w-4" />
            Import
          </Button>
          <Button
            variant="outline"
            className="hidden rounded-full border-white/15 bg-white/5 text-white/80 hover:bg-white/10 sm:flex"
          >
            <Database className="h-4 w-4" />
            Export
          </Button>
          <Button className="rounded-full bg-amber-300 font-medium text-slate-900 shadow-lg shadow-amber-400/30 hover:bg-amber-200">
            <Share2 className="h-4 w-4" />
            Share
          </Button>
        </div>
      </div>
    </header>
  )
}
