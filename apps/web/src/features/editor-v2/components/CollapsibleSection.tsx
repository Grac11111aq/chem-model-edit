import { useState } from 'react'
import { ChevronDown, ChevronRight } from 'lucide-react'

import type { ReactNode } from 'react'

interface CollapsibleSectionProps {
  title: string
  defaultOpen?: boolean
  children: ReactNode
  className?: string
  contentClassName?: string
}

export function CollapsibleSection({
  title,
  defaultOpen = false,
  children,
  className,
  contentClassName,
}: CollapsibleSectionProps) {
  const [isOpen, setIsOpen] = useState(defaultOpen)

  return (
    <div
      className={`overflow-hidden rounded-md border border-border bg-card transition-all duration-200 ${className ?? ''}`}
    >
      <button
        type="button"
        onClick={() => setIsOpen((prev) => !prev)}
        className={`flex w-full items-center gap-2 px-3 py-2 text-sm font-medium text-foreground transition-colors hover:bg-accent ${
          isOpen ? 'bg-secondary/50' : 'bg-transparent'
        }`}
      >
        {isOpen ? (
          <ChevronDown className="h-4 w-4 text-muted-foreground" />
        ) : (
          <ChevronRight className="h-4 w-4 text-muted-foreground" />
        )}
        {title}
      </button>

      {isOpen ? (
        <div
          className={`animate-in slide-in-from-top-1 border-t border-border p-3 duration-200 ${contentClassName ?? ''}`}
        >
          {children}
        </div>
      ) : null}
    </div>
  )
}
