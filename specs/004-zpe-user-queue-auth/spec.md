# Feature Specification: User Auth + Queue Target for ZPE Remote Compute

**Feature Branch**: `[stack/issue66-spec]`  
**Created**: 2026-01-28  
**Status**: Draft  
**Input**: Issue #66 (User-run compute with lightweight auth + per-user queue selection)

## User Scenarios & Tests (Required)

### User Story 1 - Open signup and session-based auth (Priority: P1)

Users can create accounts without admin approval and receive a session token used for all ZPE actions.
Session tokens are transported via `Authorization: Bearer <token>` and stored client-side in localStorage (MVP). CSRF is not applicable for header-based auth, but XSS hardening is required (CSP + input/output encoding + avoid inline scripts + short-lived tokens).

**Independent Test**: A new user can register, log in, and call `/auth/me` with a bearer token.

**Acceptance Scenarios**:
1. **Given** a new email/password, **When** `/auth/register` is called, **Then** a user is created and a session token is issued
2. **Given** valid credentials, **When** `/auth/login` is called, **Then** a session token is issued
3. **Given** a bearer token, **When** `/auth/me` is called, **Then** user profile is returned
4. **Given** a bearer token, **When** `/auth/logout` is called, **Then** the session is invalidated

---

### User Story 2 - Register a compute queue target (Priority: P1)

A logged-in user can register their compute server by consuming a short-lived enroll token and store a queue target under their account.

**Independent Test**: A user can create an enroll token, register a queue target, and see it in the target list.

**Acceptance Scenarios**:
1. **Given** a logged-in user, **When** `/calc/zpe/compute/enroll-tokens` is called, **Then** a short-lived token is issued
2. **Given** a valid enroll token, **When** `/calc/zpe/compute/servers/register` is called, **Then** a server_id is issued and stored under the user
3. **Given** existing targets, **When** `/calc/zpe/compute/targets` is called, **Then** the user’s targets are returned
4. **Given** a target, **When** `/calc/zpe/compute/targets/select` is called, **Then** it becomes the active queue for the user

---

### User Story 3 - Submit ZPE jobs to the user’s selected queue (Priority: P1)

ZPE job submission uses the authenticated user’s selected queue target; job status and results are only accessible to the job owner.

**Independent Test**: A user submits a job and can fetch status/results; another user is denied.

**Acceptance Scenarios**:
1. **Given** a logged-in user with an active queue target, **When** `/calc/zpe/jobs` is called, **Then** a job is enqueued to the selected queue
2. **Given** a job owner, **When** `/calc/zpe/jobs/{id}` is called, **Then** status is returned
3. **Given** a different user, **When** `/calc/zpe/jobs/{id}` is called, **Then** 403 is returned

---

### User Story 4 - Keep admin-only operations separate (Priority: P2)

Admin tokens are not required for standard users, and admin-only operations remain possible via `ZPE_ADMIN_TOKEN`.

**Independent Test**: A user can issue enroll tokens without admin secret in the UI.

**Acceptance Scenarios**:
1. **Given** no admin token, **When** `/calc/zpe/compute/enroll-tokens` is called with user auth, **Then** a token is issued
2. **Given** admin token, **When** `/calc/zpe/compute/enroll-tokens` is called, **Then** a token is issued even without user auth

## Requirements (Required)

### Functional Requirements

- **FR-001**: Provide open user registration (no invite/approval)
- **FR-002**: Authenticate via session tokens stored in Redis and sent as `Authorization: Bearer <token>`; tokens have >=128 bits entropy and are stored client-side in localStorage (MVP) with XSS hardening (CSP + input/output encoding + avoid inline scripts). If cookies are ever used, require HttpOnly + Secure + SameSite and CSRF tokens.
- **FR-003**: Sessions have a configurable TTL (default 7 days), are refreshed on use (sliding TTL), and are revoked by deleting the Redis session key on logout
- **FR-004**: Allow authenticated users to create compute enroll tokens
- **FR-005**: Register compute servers under a user and store a queue target (same Redis, queue name only)
- **FR-006**: Users can list and select their active queue target
- **FR-007**: ZPE job submission enqueues to the user’s selected queue target
- **FR-008**: Job status/result access is limited to the job owner
- **FR-009**: Admin token usage remains supported and is not required for end users
- **FR-010**: Password hashing uses Argon2id (m=19456 KiB, t=2, p=1) as the target; MVP may use PBKDF2-HMAC-SHA256 (>=210k iterations) until Argon2id is available; bcrypt (cost >= 10) only for legacy migration
- **FR-011**: Rate-limit /auth/register and /auth/login (per-account + per-IP), lock out after 5-10 consecutive failures, and cap ~100 failed attempts per account per hour with exponential backoff
- **FR-012**: When enroll tokens are issued via admin token, emit an audit log (actor, timestamp, source IP, target)

### Non-Functional Requirements

- **NFR-001**: No external auth provider required (no Clerk/Auth0)
- **NFR-002**: No email verification required for MVP
- **NFR-003**: API changes should not break existing ZPE compute worker

## Key Entities

- **User**: `user_id`, `email`, `password_hash`, `created_at`
- **Session**: `session_token`, `user_id`, `expires_at`
- **QueueTarget**: `target_id`, `user_id`, `queue_name`, `server_id`, `registered_at`, `name` (optional)
- **JobOwner**: `job_id` -> `user_id` mapping stored in Redis

## Success Criteria

- **SC-001**: A new user can register, log in, and submit ZPE jobs
- **SC-002**: Jobs are enqueued to the selected queue target and executed by the remote worker
- **SC-003**: A user cannot access another user’s job status or results
- **SC-004**: Admin token remains optional for enroll token issuance

## Out of Scope

- OAuth / SSO / email verification
- Cross-Redis queue targets
- Multi-tenant RBAC beyond owner-based access
