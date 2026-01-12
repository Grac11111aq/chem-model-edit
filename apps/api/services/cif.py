from __future__ import annotations

from io import BytesIO
import json
import os
from pathlib import Path
import shutil
import subprocess
from tempfile import TemporaryDirectory
from threading import Lock
from typing import Optional

from ase import Atoms as ASEAtoms
from ase.io import write as ase_write

_CONVERT_LOCK = Lock()
BCIF_TIMEOUT_SECONDS = 30


def atoms_to_cif(atoms: ASEAtoms, *, wrap: bool = False) -> str:
    buffer = BytesIO()
    ase_write(buffer, atoms, format="cif", wrap=wrap)
    return buffer.getvalue().decode("utf-8")


def _encoding_hints(precision: int) -> list[dict[str, object]]:
    columns = [
        "Cartn_x",
        "Cartn_y",
        "Cartn_z",
        "fract_x",
        "fract_y",
        "fract_z",
    ]
    return [
        {
            "categoryName": "atom_site",
            "columnName": column,
            "encoding": "delta",
            "precision": precision,
        }
        for column in columns
    ]


def _molstar_cif2bcif_path() -> Path:
    env_path = os.environ.get("MOLSTAR_CIF2BCIF_PATH")
    if env_path:
        return Path(env_path).expanduser()
    root = Path(__file__).resolve().parents[3]
    return (
        root
        / "apps"
        / "web"
        / "node_modules"
        / "molstar"
        / "lib"
        / "commonjs"
        / "cli"
        / "cif2bcif"
        / "index.js"
    )


def _node_bin() -> str:
    node_bin = os.environ.get("NODE_BIN", "node")
    resolved = shutil.which(node_bin)
    if resolved:
        return resolved
    if Path(node_bin).exists():
        return node_bin
    raise ValueError(
        "Node executable not found. Set NODE_BIN or ensure node is on PATH."
    )


def cif_to_bcif(
    cif_text: str,
    *,
    lossy: bool,
    precision: Optional[int],
) -> bytes:
    cif2bcif = _molstar_cif2bcif_path()
    if not cif2bcif.exists():
        raise ValueError(
            "molstar cif2bcif was not found. Set MOLSTAR_CIF2BCIF_PATH or install dependencies in apps/web."
        )

    with TemporaryDirectory() as tmpdir:
        tmp_path = Path(tmpdir)
        src_path = tmp_path / "input.cif"
        out_path = tmp_path / "output.bcif"
        src_path.write_text(cif_text, encoding="utf-8")

        node_bin = _node_bin()
        cmd = [node_bin, str(cif2bcif), str(src_path), str(out_path)]
        if lossy:
            if precision is None:
                precision = 3
            config_path = tmp_path / "config.json"
            config_path.write_text(
                json.dumps(_encoding_hints(precision)),
                encoding="utf-8",
            )
            cmd.extend(["-c", str(config_path)])

        with _CONVERT_LOCK:
            try:
                subprocess.run(
                    cmd,
                    check=True,
                    capture_output=True,
                    text=True,
                    timeout=BCIF_TIMEOUT_SECONDS,
                )
            except subprocess.TimeoutExpired as exc:
                raise ValueError("CIF->bcif conversion timed out.") from exc
            except subprocess.CalledProcessError as exc:
                detail = exc.stderr.strip() or exc.stdout.strip()
                raise ValueError(f"CIF->bcif conversion failed: {detail}") from exc

        if not out_path.exists():
            raise ValueError("CIF->bcif output was not found.")

        return out_path.read_bytes()
