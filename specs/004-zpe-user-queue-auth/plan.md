# Implementation Plan: User Auth + Queue Target for ZPE Remote Compute

**Branch**: `[stack/issue66-spec]` | **Date**: 2026-01-28 | **Spec**: `specs/004-zpe-user-queue-auth/spec.md`
**Input**: Issue #66 and spec above

## Summary

Introduce lightweight user authentication (open signup + Redis-backed sessions), per-user compute queue targets, and job ownership checks. ZPE job submission will use the selected queue target for the authenticated user. Admin-only enroll token issuance remains supported, but user-auth can also issue tokens for self-managed compute registration.

## Technical Context

- **Backend**: FastAPI + Redis (existing), RQ for compute queue
- **Frontend**: TanStack Start + existing ZPE UI
- **Storage**: Redis for sessions, user profiles, queue targets, and job ownership
- **Session TTL**: 7 days (sliding)
- **Password Hashing**: Argon2id (m=19456 KiB, t=2, p=1) target; MVP may use PBKDF2-HMAC-SHA256 (>=210k) until Argon2id is available. Bcrypt (cost >= 10) only for legacy migration.
- **Auth Hardening**: rate-limit /auth/register and /auth/login; lock out after 5-10 consecutive failures and cap ~100 failed attempts per account per hour (with backoff).
- **Token Transport**: Authorization Bearer header (localStorage in web client), so CSRF is not applicable; XSS hardening required.

## Project Structure

```plaintext
specs/004-zpe-user-queue-auth/
├── spec.md
├── plan.md
└── tasks.md

apps/api/
├── services/auth/         # new auth store + helpers
├── services/zpe/           # extend enroll + queue target store
└── main.py                 # add auth + queue endpoints + guards

apps/web/
├── src/lib/api.ts          # add auth + queue API
├── src/lib/auth.ts         # session storage helper
└── src/features/editor-v2/ # add login + compute settings UI
```

## Work Phases (Stacked PRs)

1. **Spec & Plan (PR-1)**
   - Add spec/plan/tasks for Issue #66

2. **Backend Auth Core (PR-2)**
   - Add user + session storage (Redis)
   - Add `/auth/register`, `/auth/login`, `/auth/logout`, `/auth/me`
   - Add session dependency + sliding TTL
   - Add rate limiting/lockout for auth endpoints (per-account + per-IP)

3. **Queue Target + ZPE Ownership (PR-3)**
   - Store queue targets per user, select active target
   - Allow user-auth to issue enroll tokens
   - Log admin-token enroll issuance (actor/time/IP/target)
   - Store job owner mapping and enforce on status/result/files
   - Enqueue jobs to selected queue target

4. **Frontend UI + Docs + Tests (PR-4)**
   - Login/register UI and session storage
   - Queue target panel (register/list/select)
   - Attach bearer token to ZPE calls
   - Add tests + docs for worker registration flow

## Rollback Plan

- Disable auth by bypassing auth dependency (all endpoints behave as anonymous)
- Invalidate existing sessions before bypassing (delete `auth:session:*` keys) or set a short global TTL; after bypass, ignore Authorization headers so stale tokens cannot access the default queue
- Fall back to single default queue in settings
- Remove owner checks from job endpoints
- Keep Redis data; no migration rollback needed
