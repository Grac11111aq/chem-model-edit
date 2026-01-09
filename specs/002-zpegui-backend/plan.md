# Implementation Plan: ZPEGUI バックエンド移植

**Branch**: `[dev]` | **Date**: 2026-01-07 | **Spec**: `specs/002-zpegui-backend/spec.md`
**Input**: Feature specification from `specs/002-zpegui-backend/spec.md`

## Summary

QE .in 解析と ZPE/Vibrations 計算を FastAPI + RQ で非同期実行できるようにする。  
API は `/calc/zpe/parse` と `/calc/zpe/jobs` 系を追加し、ワーカーで ASE Vibrations + QE を実行。結果は JSON + CSV/summary で返却する。

## Technical Context

- **Language/Version**: Python 3.13
- **Primary Dependencies**: FastAPI, ASE, pymatgen, RQ, redis-py, pydantic-settings
- **Storage**: Redis（ジョブ状態） + ファイル（作業ディレクトリ/結果）
- **Testing**: pytest（可能なら fakeredis で job キューを検証）
- **Target Platform**: Linux server
- **Project Type**: API backend only
- **Performance Goals**: parse 1k atoms < 2s / job enqueue < 500ms
- **Constraints**: `pw.x` / pseudopotentials がサーバに存在する前提

## Project Structure

### Documentation (this feature)

```
specs/002-zpegui-backend/
├── spec.md
├── plan.md
└── tasks.md
```

### Source Code (repository root)

```
apps/api/
├── main.py              # API routes (add /calc/zpe/*)
├── models.py            # request/response models
├── services/
│   ├── zpe.py           # 解析/計算ユーティリティ
│   └── zpe_worker.py    # RQジョブ実体
├── worker.py            # RQ worker entrypoint
└── tests/
    ├── test_zpe_parse.py
    └── test_zpe_jobs.py
```

**Structure Decision**: 既存の `apps/api` に ZPE 機能を追加し、RQワーカーは同一プロジェクト内に配置する。

## Work Phases

1. **Design/Setup**: 依存追加、設定（Redis URL / 仕事ディレクトリ / QE実行 / pseudo_dir）
2. **Core Implementation**: 解析/ジョブ投入/ワーカー/結果取得の実装
3. **Verification**: 単体テスト追加、pytest/mypy/ruff 実行

## Rollback Plan

- 新規 API ルートを削除し、`services/zpe*.py` と worker を撤去
- 依存追加を `pyproject.toml` から戻す
- 既存 API には影響を与えないため影響範囲は限定的
