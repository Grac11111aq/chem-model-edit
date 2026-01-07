import { createFileRoute } from '@tanstack/react-router'

import EditorPage from '@/features/editor/EditorPage'

export const Route = createFileRoute('/editor')({
  component: EditorPage,
})
