declare module 'molstar/lib/apps/viewer/app' {
  export class Viewer {
    constructor(container: HTMLElement, options?: Record<string, unknown>)
    static create(
      container: HTMLElement,
      options?: Record<string, unknown>,
    ): Promise<Viewer>
    loadStructureFromData(data: string, format: string): void
    dispose(): void
    plugin: any
  }
}

declare module 'molstar/lib/mol-plugin-ui/skin/light.scss'
