from __future__ import annotations

import os
import sys
from pathlib import Path


def can_rename(base: Path) -> bool:
    try:
        base.mkdir(parents=True, exist_ok=True)
        src_dir = base / ".uv-rename-a"
        dst_dir = base / ".uv-rename-b"
        src_dir.mkdir(parents=True, exist_ok=True)
        dst_dir.mkdir(parents=True, exist_ok=True)
        src = src_dir / "file.txt"
        src.write_text("x", encoding="utf-8")
        dest = dst_dir / "file.txt"
        os.rename(src, dest)
        for path in (src_dir, dst_dir):
            for child in path.glob("*"):
                child.unlink(missing_ok=True)
            path.rmdir()
        return True
    except Exception:
        return False


def main() -> int:
    candidates = [c for c in sys.argv[1:] if c]
    for candidate in candidates:
        path = Path(os.path.expanduser(candidate))
        if can_rename(path):
            print(path)
            return 0
    if candidates:
        print(candidates[0])
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
