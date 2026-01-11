/* @vitest-environment jsdom */
import { cleanup, fireEvent, render, screen } from '@testing-library/react'
import { afterEach, describe, expect, it } from 'vitest'

import TransplantPage from './TransplantPage'

afterEach(() => {
  cleanup()
})

describe('TransplantPage', () => {
  it('renders heading and action button', () => {
    render(<TransplantPage />)

    expect(screen.getByText('Delta Transplant')).toBeTruthy()
    expect(screen.getByRole('button', { name: /run transplant/i })).toBeTruthy()
  })

  it('shows an error when required inputs are missing', () => {
    render(<TransplantPage />)

    fireEvent.click(screen.getByRole('button', { name: /run transplant/i }))
    expect(
      screen.getByText(
        'small .in / small .out / large .in の3つを入力してください。',
      ),
    ).toBeTruthy()
  })
})
