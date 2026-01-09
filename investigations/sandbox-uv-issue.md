# Codexサンドボックスでのuv/just実行問題

## 問題
- Codex CLI の `sandbox_mode=workspace-write` 環境で `uv sync` が失敗する。
- `just api` / `just dev` 実行時に uv が内部で行う rename が `Invalid cross-device link (os error 18)` になり停止。

## 原因
- サンドボックス内でディレクトリ跨ぎの rename が制限されている。
- `writable_roots` や `--add-dir` は「書き込み許可」の追加であり、rename 制限（EXDEV）は解消しない。

## やったこと
- 新規 worktree を作成して再現確認。
- `uv sync` 実行時に `.uv-cache/.tmp* -> .uv-cache/archive-v0/*` の rename で失敗することを確認。
- `--no-cache` でも `.uv-tmp/.tmp* -> .uv-tmp/*/archive-v0/*` の rename で失敗することを確認。
- `os.rename` の簡易テストで、同一ディレクトリ内の rename は成功・ディレクトリ跨ぎで EXDEV になることを確認。
- `UV_CACHE_DIR` / `TMPDIR` を `/home/grace/.codex` 配下に固定し、`writable_roots` / `--add-dir` を試すも EXDEV は解消せず。
- `danger-full-access` または昇格（on-request）では `uv sync` が成功することを確認。

## 結果
- 根本原因はサンドボックスの rename 制限（EXDEV）。
- `writable_roots` 追加のみでは uv の rename 問題を解決できない。

## 対応方針
- uv を使う作業は昇格（on-request）で実行する。
- もしくは `danger-full-access` で作業する。
- `Justfile` は `/home/grace/.codex/uv-cache` / `/home/grace/.codex/uv-tmp` を参照する前提で運用（ただし昇格が必要）。
