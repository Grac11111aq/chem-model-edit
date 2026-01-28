import type { AuthSession } from './types'

const SESSION_KEY = 'chem-model-auth-session'

export function getStoredSession(): AuthSession | null {
  if (typeof window === 'undefined') {
    return null
  }
  const raw = window.localStorage.getItem(SESSION_KEY)
  if (!raw) {
    return null
  }
  try {
    return JSON.parse(raw) as AuthSession
  } catch (_err) {
    return null
  }
}

export function storeSession(session: AuthSession) {
  if (typeof window === 'undefined') {
    return
  }
  window.localStorage.setItem(SESSION_KEY, JSON.stringify(session))
}

export function clearSession() {
  if (typeof window === 'undefined') {
    return
  }
  window.localStorage.removeItem(SESSION_KEY)
}

export function getAuthToken(): string | null {
  const session = getStoredSession()
  return session?.token ?? null
}
