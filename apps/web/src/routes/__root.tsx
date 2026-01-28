import { Suspense, lazy } from 'react'
import { Outlet, createRootRoute } from '@tanstack/react-router'

import Header from '../components/Header'

const RouterDevtoolsPanel = import.meta.env.DEV
  ? lazy(() =>
      import('@tanstack/react-router-devtools').then((module) => ({
        default: module.TanStackRouterDevtoolsPanel,
      })),
    )
  : null

const Devtools = import.meta.env.DEV
  ? lazy(() =>
      import('@tanstack/react-devtools').then((module) => ({
        default: module.TanStackDevtools,
      })),
    )
  : null

export const Route = createRootRoute({
  component: RootLayout,
})

function RootLayout() {
  return (
    <>
      <Header />
      <Outlet />
      {Devtools && RouterDevtoolsPanel ? (
        <Suspense fallback={null}>
          <Devtools
            config={{
              position: 'bottom-right',
            }}
            plugins={[
              {
                name: 'Tanstack Router',
                render: <RouterDevtoolsPanel />,
              },
            ]}
          />
        </Suspense>
      ) : null}
    </>
  )
}
