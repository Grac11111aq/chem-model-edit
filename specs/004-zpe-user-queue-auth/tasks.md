# Task List: User Auth + Queue Target for ZPE Remote Compute

**Input**: `specs/004-zpe-user-queue-auth/spec.md` / `specs/004-zpe-user-queue-auth/plan.md`

## PR-1: Spec/Plan
- [ ] T401 [P] Add spec/plan/tasks for Issue #66

## PR-2: Backend Auth Core
- [ ] T410 [P] Implement user + session store in Redis
- [ ] T411 [P] Add /auth/register, /auth/login, /auth/logout, /auth/me
- [ ] T412 [P] Add auth dependency with sliding TTL

## PR-3: Queue Target + Ownership
- [ ] T420 [P] Add queue target store and selection endpoints
- [ ] T421 [P] Allow user-auth to issue enroll tokens
- [ ] T422 [P] Store job owner mapping and enforce on ZPE endpoints
- [ ] T423 [P] Enqueue ZPE jobs to user-selected queue target

## PR-4: Frontend + Docs + Tests
- [ ] T430 [P] Add login/register UI + session storage
- [ ] T431 [P] Add queue target UI for register/list/select
- [ ] T432 [P] Attach auth token to ZPE API calls
- [ ] T433 [P] Add docs for worker registration flow
- [ ] T434 [P] Add tests for auth + queue target + owner checks
