import { devtools } from '@tanstack/devtools-vite'
import { tanstackRouter } from '@tanstack/router-plugin/vite'
import viteReact from '@vitejs/plugin-react'
import viteTsConfigPaths from 'vite-tsconfig-paths'
import tailwindcss from '@tailwindcss/vite'
import { defineConfig } from 'vite'
import type { PluginOption } from 'vite'

const config = defineConfig(({ mode }) => {
  const isTest = mode === 'test' || process.env.VITEST === 'true'
  const plugins = [
    !isTest &&
      devtools({
        eventBusConfig: {
          enabled: process.env.TANSTACK_DEVTOOLS_EVENTBUS === 'true',
          port: Number(process.env.TANSTACK_DEVTOOLS_PORT ?? 42069),
        },
      }),
    !isTest && tanstackRouter({ target: 'react' }),
    // this is the plugin that enables path aliases
    viteTsConfigPaths({
      projects: ['./tsconfig.json'],
    }),
    tailwindcss(),
    viteReact(),
  ].filter(Boolean) as Array<PluginOption>

  return {
    plugins,
    optimizeDeps: {
      include: ['dockview', 'dockview-react', 'dockview-core'],
    },
  }
})

export default config
