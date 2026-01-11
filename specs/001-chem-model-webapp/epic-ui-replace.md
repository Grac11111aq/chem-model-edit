# Epic: Editor UI Replace (editor-v2)

## Goal
Build a new editor experience at `/editor-v2` and iterate toward the final UI, without blocking ongoing development. The new UI should center on **immutable structures**: imported `.in` data is treated as read-only, and tools generate new structures instead of in-place editing.

## Principles
- **Immutable inputs**: imported files are preserved; no direct edits to original data.
- **Generated outputs**: tools (Transfer / Supercell / etc.) create new structures.
- **Provenance-first**: each structure should record how it was created (source + operation).
- **Low-risk rollout**: existing `/editor` can remain unchanged; parity is not required.

## Scope
- New UI layout and navigation at `/editor-v2`.
- File Manager lists structures + creation lineage.
- Workspace panels render per-structure views.
- Tool panels provide structure-generating actions.
- ZPE/Vibrations stays **preview-only** for now.

## Non-goals (for this epic)
- Full parity with `/editor`.
- Rewriting backend APIs.
- ZPE/Vibrations feature wiring.

## Definition of Done
- `/editor-v2` is usable and matches the new mental model.
- Structure creation flows are centered on **generation** instead of direct edits.
- Existing `/editor` does not need full functional parity.

## Progress Log
- 2026-01-10: `/editor-v2` route + mock layout + tool renaming + vibration preview shell.

## Notes
Feature Issues will be created only as needed during UI review and iteration.
