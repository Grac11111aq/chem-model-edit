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

describe('SupercellPage', () => {
  it('renders headline and generate button', async () => {
    const { default: SupercellPage } = await import('./SupercellPage')
    render(<SupercellPage />)

    expect(screen.getByText('Supercell Builder')).toBeTruthy()
    expect(screen.getByRole('button', { name: /generate/i })).toBeTruthy()
  }, 15000)

  it('shows an error for invalid tile pattern', async () => {
    const { default: SupercellPage } = await import('./SupercellPage')
    render(<SupercellPage />)

    fireEvent.change(screen.getByLabelText('Tile pattern'), {
      target: { value: 'AC' },
    })
    expect(
      screen.getByText('タイルは A/B のみで指定してください。'),
    ).toBeTruthy()
  })
})
