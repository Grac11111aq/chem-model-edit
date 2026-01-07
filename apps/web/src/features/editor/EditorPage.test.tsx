/* @vitest-environment jsdom */
import { cleanup, fireEvent, render, screen } from '@testing-library/react'
import { afterEach, describe, expect, it, vi } from 'vitest'

vi.mock('@/components/molstar/MolstarViewer', () => ({
  default: () => <div data-testid="molstar-viewer" />,
}))

const globalWithResizeObserver = globalThis as unknown as {
  ResizeObserver?: typeof ResizeObserver
}

if (!globalWithResizeObserver.ResizeObserver) {
  globalWithResizeObserver.ResizeObserver = class {
    observe() {}
    unobserve() {}
    disconnect() {}
  }
}

afterEach(() => {
  cleanup()
})

describe('EditorPage', () => {
  it('renders key sections', async () => {
    const { default: EditorPage } = await import('./EditorPage')
    render(<EditorPage />)

    expect(screen.getByText('Atom Table')).toBeTruthy()
    expect(screen.getByText('Mol* Preview')).toBeTruthy()
    expect(screen.getByText('Structures')).toBeTruthy()
  }, 15000)

  it('adds an atom row', async () => {
    const { default: EditorPage } = await import('./EditorPage')
    render(<EditorPage />)

    fireEvent.click(screen.getByRole('button', { name: /add atom/i }))
    expect(screen.getByLabelText('atom 1 symbol')).toBeTruthy()
  })
})
