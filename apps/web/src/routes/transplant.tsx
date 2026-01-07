import { createFileRoute } from '@tanstack/react-router'

import TransplantPage from '@/features/transplant/TransplantPage'

export const Route = createFileRoute('/transplant')({
  component: TransplantPage,
})
