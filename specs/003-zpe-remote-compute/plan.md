# Implementation Plan: ZPE 計算サーバ分離

**Branch**: `[refactor/zpe-remote-compute]` | **Date**: 2026-01-10 | **Spec**: `specs/003-zpe-remote-compute/spec.md`
**Input**: Feature specification from `specs/003-zpe-remote-compute/spec.md`

## Summary

ZPE 計算を `control-plane`（API）と `compute-plane`（計算）に分離し、サービス提供側では計算が走らない構成にする。ジョブは Redis/RQ で計算サーバーに委譲し、結果は共有ストア（初期は Redis）で受け渡す。E2E/CI 向けにモック計算モードも提供する。

## Technical Context

- **Language/Version**: Python 3.13
- **Primary Dependencies**: FastAPI, RQ, redis-py, pydantic-settings
- **Result Store**: Public Redis + TLS + ACL（将来 S3 へ差し替え可能な設計）
- **Testing**: pytest（モックバックエンド + ストアのユニットテスト）
- **Target Platform**: API: 任意のクラウド / Compute: 計算サーバー
- **Constraints**: API 側は `pw.x` や pseudo を持たない前提

## Project Structure

### Documentation (this feature)

```
specs/003-zpe-remote-compute/
├── spec.md
├── plan.md
└── tasks.md
```

### Source Code (repository root)

```
apps/api/
├── services/zpe/
│   ├── backends.py      # 計算委譲アダプタ
│   ├── result_store.py  # 結果保存/取得
│   └── settings.py      # compute mode など追加
└── main.py              # /calc/zpe/* でバックエンド切替を利用
```

**Structure Decision**: 既存の `apps/api/services/zpe` にアダプタ層を追加し、API 自体は `/calc/zpe/*` を維持する。

## Work Phases

1. **Design/Setup**
   - `ZPE_COMPUTE_MODE` / `ZPE_RESULT_STORE` など設定追加
   - Public Redis + TLS + ACL 前提の構成を明示
   - 短期登録トークン方式（将来のユーザー単位登録へ拡張）を設計
   - 結果ストアとバックエンドのインターフェース定義
2. **Core Implementation**
   - `remote-queue`（RQ + Redis）バックエンド
   - `mock` バックエンド（決定論的なダミー結果）
   - 結果保存/取得の統一（Redis）
3. **Verification**
   - 共有ストア経由の結果取得テスト
   - モックモードのE2E向け検証

## Rollback Plan

- 追加したバックエンド/ストア層を削除
- `ZPE_COMPUTE_MODE` を `local` に戻す
- `/calc/zpe/*` の挙動を従来のローカル方式に戻す
