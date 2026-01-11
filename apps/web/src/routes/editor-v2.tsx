import { createFileRoute } from '@tanstack/react-router'

import EditorV2Page from '@/features/editor-v2/EditorV2Page'

export const Route = createFileRoute('/editor-v2')({
  component: EditorV2Page,
})
