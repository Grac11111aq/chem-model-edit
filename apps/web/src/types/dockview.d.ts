import 'dockview-react'

declare module 'dockview-react' {
  interface DockviewApi {
    activePanel?: unknown
    onDidAddPanel: (listener: (panel: unknown) => void) => {
      dispose: () => void
    }
    onDidRemovePanel: (listener: (panel: unknown) => void) => {
      dispose: () => void
    }
  }
}
