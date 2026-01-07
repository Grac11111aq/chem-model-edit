import { createFileRoute } from '@tanstack/react-router'

import SupercellPage from '@/features/supercell/SupercellPage'

export const Route = createFileRoute('/supercell')({
  component: SupercellPage,
})
