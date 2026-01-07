# -*- coding: utf-8 -*-
"""
ZPE GUI (QE + ASE Vibrations)
- QEインプットを読み込んで、表＋3Dで可動原子(indices)を手動選択
- [計算]ボタンで選択原子のみ有限差分振動 → 周波数・ZPE算出
- 結果をテキスト/CSVに保存（タイムスタンプ付き）
- EspressoProfile のAPI差を自動吸収（argv/command/pseudo_dir など）
- 例外はGUIを固めずに安全に表示（lambda のデフォルト引数で NameError を回避）
- 計算前に Vibrations キャッシュ(vib_name.*)を“ファイル単位”で健全性チェックし、0バイト/壊れた JSON だけ削除してから再開
- 新機能:
   * 計算開始時間と終了時間を出力ファイルに記録
   * Environ使用の有無をチェックボックスで指定し出力に記録
   * [新規計算]ボタンで既存計算結果を削除して新規計算か継続計算かを出力に記録
"""

import os, re, shutil, inspect, threading, csv, traceback, json
from glob import glob
from pathlib import Path
from datetime import datetime

import numpy as np
import tkinter as tk
from tkinter import ttk, messagebox, filedialog

from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
import matplotlib
matplotlib.use("TkAgg")
# 日本語フォント設定
matplotlib.rcParams['font.family'] = 'Noto Sans CJK JP'   # インストールしたフォント名

import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

from ase.io import read
from ase.constraints import FixAtoms
from ase.calculators.espresso import Espresso, EspressoProfile
from ase.vibrations import Vibrations

from ase.units import kB   # ASE からはボルツマン定数だけ使う
h = 6.62607015e-34         # J*s  (2023 CODATA 国際定義)

QE_IN_DEFAULT = "WO3_slab_202_Vo-OH_8A_scf_k311.in"
PSEUDOS = {"W": "W.upf", "O": "O.upf", "H": "H.upf"}
OUTDIR = "calc_scratch"
VIB_NAME = "vib_ads"
DELTA = 0.01
KPTS = (3, 3, 1)
LOW_CUT_CM = 50.0
USE_MPI = True
NP_CORE = 12

import math

# --- QE namelist & K_POINTS(automatic) パーサ ---
def _qe_parse_namelist(qe_path: str, name: str) -> dict:
    """QE入力から &name ... / のブロックを読み取り dict で返す。存在しなければ {}。
    - コメント(!以降)除去、値は Python 型へ: .true./.false. → bool, 数値 → int/float, 文字列 → str
    """
    txt = Path(qe_path).read_text(encoding="utf-8", errors="ignore")
    m = re.search(rf"&\s*{name}\b(.*?)\n\s*/", txt, re.I | re.S)
    if not m:
        return {}
    body = m.group(1)
    out = {}
    for line in body.splitlines():
        line = line.split('!')[0].strip().rstrip(',')
        if not line or '=' not in line:
            continue
        k, v = line.split('=', 1)
        key = k.strip().lower()
        val = v.strip()
        lval = val.lower()
        if lval in ('.true.', 'true'):
            py = True
        elif lval in ('.false.', 'false'):
            py = False
        else:
            try:
                if any(ch in lval for ch in ['.', 'e', 'd']):
                    py = float(lval.replace('d', 'e'))
                    if math.isfinite(py) and abs(py - int(py)) < 1e-15:
                        py = int(py)
                else:
                    py = int(lval)
            except Exception:
                py = val.strip().strip('"\'')
        out[key] = py
    return out

def _qe_parse_kpoints_automatic(qe_path: str):
    """K_POINTS automatic を ((nkx,nky,nkz), (sx,sy,sz)) で返す。無ければ None。"""
    txt = Path(qe_path).read_text(encoding="utf-8", errors="ignore")
    m = re.search(
        r"K_POINTS\s*\(\s*automatic\s*\)\s*\n\s*([-\d]+)\s+([-\d]+)\s+([-\d]+)\s+([-\d]+)\s+([-\d]+)\s+([-\d]+)",
        txt, re.I)
    if not m:
        return None
    nkx, nky, nkz, sx, sy, sz = (int(m.group(i)) for i in range(1, 7))
    return (nkx, nky, nkz), (sx, sy, sz)

def _prepare_environ_files(vib_name: str, n_mobile: int, environ_src: Path):
    """Vibrations が使う全ディレクトリに environ.in を前置きする。
    対象: vib_name.eq と vib_name.{i}{axis}{±} for i=0..n_mobile-1, axis∈{x,y,z}
    """
    # eq
    eq_dir = Path(f"{vib_name}.eq")
    eq_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy2(environ_src, eq_dir / "environ.in")

    # 全変位（+/- × xyz × n_mobile）
    for i in range(n_mobile):
        for axis in ("x", "y", "z"):
            for sign in ("+", "-"):
                d = Path(f"{vib_name}.{i}{axis}{sign}")
                d.mkdir(parents=True, exist_ok=True)
                shutil.copy2(environ_src, d / "environ.in")
    # 追加：calc.directory（= "vib_ads_vac"）配下にも同じ階層ができることがあるので、そこにも配る
    alt_root = Path("vib_ads_vac")
    if alt_root.exists():
        # eq
        (alt_root / f"{vib_name}.eq").mkdir(parents=True, exist_ok=True)
        shutil.copy2(environ_src, alt_root / f"{vib_name}.eq" / "environ.in")
        # 全変位
        for i in range(n_mobile):
            for axis in ("x", "y", "z"):
                for sign in ("+", "-"):
                    d = alt_root / f"{vib_name}.{i}{axis}{sign}"
                    d.mkdir(parents=True, exist_ok=True)
                    shutil.copy2(environ_src, d / "environ.in")

# Vibrations のキャッシュを“ファイル単位”で健全性チェックする関数
def sanitize_vib_cache(name: str, natoms: int | None = None) -> dict:
    """
    Vibrations キャッシュの健全性チェックを行い、壊れた JSON を個別削除する。
    対象:
      - ./{name}.eq, ./{name}.*
      - ./vib_ads_vac/{name}.eq, ./vib_ads_vac/{name}.*
      - ./{name}/ 以下（※ここは「非ゼロなら削除しない」= 緩モード）
      - ./vib_ads_vac/{name}/ 以下（※ここは厳格モード）

    ルール（厳格モード = 通常）:
      - 0バイト / JSON失敗 / 'forces' 欠損・None・空 → 削除
      - natoms 指定時に forces 長さ不一致 → 削除
      - forces に NaN/inf → 削除
      - forces 形状が (*,3) でない → 削除

    ルール（緩モード = {name}/ 以下のみ）:
      - 0バイト / JSON失敗 / 'forces' 欠損・None・空 → 削除
      - 上記以外（長さ不一致や NaN/形状不一致）は「残す」
        → ご要望通り「0KBでない json は削除しない」方針
    """
    roots_strict = [
        Path("."),                  # 直下の vib_ads.* / vib_ads.eq
        Path("vib_ads_vac"),       # vib_ads_vac の配下
        Path("vib_ads_vac") / name # vib_ads_vac/vib_ads の再帰
    ]
    roots_relaxed = [
        Path(name)                  # vib_ads/ 配下（再帰）：ここは緩モード
    ]

    deleted = 0
    checked = 0
    deleted_list = []

    def check_and_maybe_delete(jpath: Path, strict: bool):
        nonlocal deleted, checked, deleted_list
        checked += 1
        try:
            if not jpath.exists():
                return
            # 0バイトは無条件削除（共通）
            if jpath.stat().st_size == 0:
                raise ValueError("empty cache")
            # JSON を開く
            with jpath.open("r", encoding="utf-8") as f:
                data = json.load(f)
            forces = data.get("forces", None)

            # 'forces' 欠損/空/None は削除（共通）
            if forces is None or (isinstance(forces, list) and len(forces) == 0):
                raise ValueError("forces missing/empty")

            # ここからは厳格モードでのみ判定
            if strict:
                if natoms is not None and len(forces) != natoms:
                    raise ValueError(f"forces length {len(forces)} != natoms {natoms}")
                arr = np.array(forces, dtype=float)
                if arr.ndim != 2 or arr.shape[1] != 3:
                    raise ValueError(f"forces bad shape {arr.shape}")
                if not np.isfinite(arr).all():
                    raise ValueError("forces contains NaN/inf")

        except Exception:
            try:
                jpath.unlink()
                deleted += 1
                deleted_list.append(str(jpath))
            except Exception:
                pass  # 権限などで消せない場合はスキップ

    # 1) eq / 変位ディレクトリ（非再帰）：厳格
    for root in (Path("."), Path("vib_ads_vac")):
        eq_json = root / f"{name}.eq" / "cache.eq.json"
        if eq_json.exists():
            check_and_maybe_delete(eq_json, strict=True)
        for d in root.glob(f"{name}.*"):
            if d.is_dir():
                for jf in d.glob("cache*.json"):
                    check_and_maybe_delete(jf, strict=True)

    # 2) 再帰スキャン：vib_ads_vac/{name}/ は厳格、{name}/ は緩
    #   - vib_ads/ 配下（緩モード）：0KB/JSON失敗/forces欠損のみ消す
    for base in roots_relaxed:
        if base.exists() and base.is_dir():
            for jf in base.rglob("cache*.json"):
                check_and_maybe_delete(jf, strict=False)

    #   - vib_ads_vac/vib_ads/ 配下（厳格モード）
    for base in (Path("vib_ads_vac") / name,):
        if base.exists() and base.is_dir():
            for jf in base.rglob("cache*.json"):
                check_and_maybe_delete(jf, strict=True)

    return {"deleted_files": deleted, "checked_files": checked, "deleted_list": deleted_list}

def read_qe_atoms(qe_in: str):
    atoms = read(qe_in, format='espresso-in')
    with open(qe_in, "r", encoding="utf-8") as f:
        txt = f.read()
    m = re.search(r"ATOMIC_POSITIONS[^\n]*\n(.*?)\n\s*CELL_PARAMETERS", txt, re.S)
    if not m:
        raise RuntimeError("ATOMIC_POSITIONS block not found in QE input.")
    apos_block = m.group(1).strip().splitlines()
    fixed_idx = []
    for i, line in enumerate(apos_block):
        if re.search(r"\s0\s+0\s+0\s*$", line):
            fixed_idx.append(i)
    if fixed_idx:
        atoms.set_constraint(FixAtoms(indices=fixed_idx))
    return atoms, fixed_idx


def build_espresso_profile(pseudo_dir: str, use_mpi: bool, np_core: int, environ: bool):
    pw = shutil.which("pw.x")
    if pw is None:
        raise SystemExit("ERROR: 'pw.x' not found.")
    argv = []
    if use_mpi:
        argv += ["mpirun", "-np", str(np_core), pw]
    else:
        argv += [pw]
    if environ:
        argv += ["-environ"]
    sig = inspect.signature(EspressoProfile)
    params = set(sig.parameters)
    if {"command", "pseudo_dir", "argv"} <= params:
        return EspressoProfile(command=("mpirun" if use_mpi else pw), pseudo_dir=pseudo_dir, argv=(argv[2:] if use_mpi else argv[1:]))
    elif {"command", "pseudo_dir"} <= params:
        return EspressoProfile(command=" ".join(argv), pseudo_dir=pseudo_dir)
    elif {"argv"} <= params:
        return EspressoProfile(argv=argv)
    else:
        raise SystemExit("Unsupported EspressoProfile signature")


class ZPEGUI(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("ZPE: select mobile atoms (Table + 3D)")
        self.geometry("1150x760")

        top = ttk.Frame(self)
        top.pack(side=tk.TOP, fill=tk.X, padx=8, pady=6)

        ttk.Label(top, text="QE input:").pack(side=tk.LEFT)
        self.path_var = tk.StringVar(value=QE_IN_DEFAULT)
        self.entry = ttk.Entry(top, width=60, textvariable=self.path_var)
        self.entry.pack(side=tk.LEFT, padx=6)
        ttk.Button(top, text="参照…", command=self.browse_file).pack(side=tk.LEFT, padx=4)
        ttk.Button(top, text="再読込", command=self.load_structure).pack(side=tk.LEFT, padx=4)

        # Environ使用チェックボックス
        self.use_environ_var = tk.BooleanVar(value=False)
        ttk.Checkbutton(top, text="Environ使用", variable=self.use_environ_var).pack(side=tk.LEFT, padx=4)

        # 新規計算ボタン
        self.btn_newcalc = ttk.Button(top, text="新規計算", command=lambda: self.run_calculation(new_calc=True), state=tk.DISABLED)
        self.btn_newcalc.pack(side=tk.RIGHT, padx=4)

        self.btn_calc = ttk.Button(top, text="継続計算", command=lambda: self.run_calculation(new_calc=False), state=tk.DISABLED)
        self.btn_calc.pack(side=tk.RIGHT, padx=4)

        left = ttk.Frame(self)
        left.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=(8,4), pady=6)
        columns = ("index","elem","x","y","z","fixed")
        self.tree = ttk.Treeview(left, columns=columns, show="headings", selectmode="extended")
        for c, w in zip(columns, (60,60,100,100,100,70)):
            self.tree.heading(c, text=c)
            self.tree.column(c, width=w, anchor="center")
        self.tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        self.tree.bind("<<TreeviewSelect>>", self.on_table_select)
        yscroll = ttk.Scrollbar(left, orient="vertical", command=self.tree.yview)
        self.tree.configure(yscroll=yscroll.set)
        yscroll.pack(side=tk.RIGHT, fill=tk.Y)

        right = ttk.Frame(self)
        right.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=(4,8), pady=6)
        self.fig = plt.Figure(figsize=(6,6))
        self.ax = self.fig.add_subplot(111, projection='3d')
        self.canvas = FigureCanvasTkAgg(self.fig, master=right)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        self.scatter = None
        self.cid_pick = None
        #ステータスバー
        bottom = ttk.Frame(self)
        bottom.pack(side=tk.BOTTOM, fill=tk.X, padx=8, pady=6)
        self.status_var = tk.StringVar(value="準備完了")
        ttk.Label(bottom, textvariable=self.status_var, anchor="w").pack(side=tk.LEFT)

        self.atoms = None
        self.fixed_idx = []
        self.selected = set()

        self.load_structure()

    def browse_file(self):
        f = filedialog.askopenfilename(filetypes=[("QE input","*.in *.pw *.qe *.txt"),("All","*.*")])
        if f:
            self.path_var.set(f)
            self.load_structure()

    def load_structure(self):
        try:
            qe_in = self.path_var.get()
            self.atoms, self.fixed_idx = read_qe_atoms(qe_in)
        except Exception as e:
            messagebox.showerror("読み込み失敗", str(e))
            return
        self.selected.clear()
        self.populate_table()
        self.draw_3d()
        self.btn_calc.config(state=tk.NORMAL)
        self.btn_newcalc.config(state=tk.NORMAL)

    def populate_table(self):
        self.tree.delete(*self.tree.get_children())
        pos = self.atoms.get_positions()
        syms = self.atoms.get_chemical_symbols()
        for i, (s, r) in enumerate(zip(syms, pos)):
            fixed = ("Yes" if i in self.fixed_idx else "No")
            self.tree.insert("", "end", iid=str(i), values=(i, s, f"{r[0]:.4f}", f"{r[1]:.4f}", f"{r[2]:.4f}", fixed))

    def draw_3d(self):
        self.ax.clear()
        pos = self.atoms.get_positions()
        colors = []
        for i in range(len(self.atoms)):
            if i in self.selected:
                colors.append("yellow")
            elif i in self.fixed_idx:
                colors.append("#c0c0c0")
            else:
                colors.append("tab:blue")
        xyz = pos
        self.scatter = self.ax.scatter(xyz[:,0], xyz[:,1], xyz[:,2], s=40, picker=True, depthshade=False, c=colors)
        self.ax.set_xlabel("x (Å)"); self.ax.set_ylabel("y (Å)"); self.ax.set_zlabel("z (Å)")
        self.ax.set_title("3D: 点クリックで選択/解除")
        if self.cid_pick is not None:
            self.canvas.mpl_disconnect(self.cid_pick)
        self.cid_pick = self.canvas.mpl_connect("pick_event", self.on_pick)
        self.canvas.draw_idle()

    def on_pick(self, event):
        if event.artist != self.scatter:
            return
        ind = event.ind
        if len(ind) == 0:   # ←ここを修正:
            return
        idx = int(ind[0])
        if idx in self.selected:
            self.selected.remove(idx)
        else:
            self.selected.add(idx)
        if str(idx) in self.tree.selection():
            self.tree.selection_remove(str(idx))
        else:
            self.tree.selection_add(str(idx)); self.tree.see(str(idx))
        self.draw_3d()

    def on_table_select(self, _evt=None):
        self.selected = set(int(i) for i in self.tree.selection())
        self.draw_3d()

    def run_calculation(self, new_calc=False):
        if self.atoms is None:
            return
        mobile = sorted(int(i) for i in self.selected if i not in self.fixed_idx)
        if not mobile:
            messagebox.showwarning("選択が空", "固定原子を除いた '可動原子' が 1 つも選ばれていません。")
            return
        if new_calc:
            for d in glob(f"{VIB_NAME}.*"):
                shutil.rmtree(d, ignore_errors=True)
        profile = build_espresso_profile(pseudo_dir=os.path.abspath('.'), use_mpi=USE_MPI, np_core=NP_CORE, environ=self.use_environ_var.get())
        self.btn_calc.config(state=tk.DISABLED)
        self.btn_newcalc.config(state=tk.DISABLED)
        start_time = datetime.now()
        #ステータスバー
        self.status_var.set("計算中…")
        threading.Thread(target=self._calc_worker, args=(mobile, profile, start_time, new_calc), daemon=True).start()

    def _calc_worker(self, mobile, profile: EspressoProfile, start_time: datetime, new_calc: bool):
        try:
            Path(OUTDIR).mkdir(parents=True, exist_ok=True)
            cleaned = sanitize_vib_cache(VIB_NAME, natoms=len(self.atoms))
            use_environ = bool(self.use_environ_var.get())
            if use_environ:
                qe_dir = Path(self.path_var.get()).resolve().parent
                cand = [qe_dir / "environ.in", Path.cwd() / "environ.in"]
                environ_src = next((p for p in cand if p.exists()), None)
                if environ_src is None:
                    raise FileNotFoundError(
                        "Environ使用がONですが environ.in が見つかりません。\n"
                         f"次のいずれかに置いてください:\n  - {qe_dir}/environ.in\n  - {Path.cwd()}/environ.in"
                    )
                # 新規/継続に関わらず、全変位ディレクトリ（eq, ±x/y/z×各モバイル）へ前置き
                _prepare_environ_files(VIB_NAME, n_mobile=len(mobile), environ_src=environ_src)

            # --- &system / &electrons を QE 入力から読み取り、ASE へ渡す ---
            qe_path = self.path_var.get()
            sys_from_in = _qe_parse_namelist(qe_path, 'system')
            ele_from_in = _qe_parse_namelist(qe_path, 'electrons')

            # 不足したときのフォールバック（必要最小限だけ補う）
            sys_defaults = {
                'occupations': 'smearing',
                'smearing': 'mp',
                'degauss': 0.02,
            }
            ele_defaults = {
                'conv_thr': 1.0e-6,
                'mixing_beta': 0.3,
                'mixing_mode': 'plain',
                'electron_maxstep': 100,
            }
            system_dict = {**sys_defaults, **sys_from_in}
            electrons_dict = {**ele_defaults, **ele_from_in}

            # K_POINTS automatic があればそれを使い、無ければ既定 KPTS
            kp = _qe_parse_kpoints_automatic(qe_path)
            kpts_use = KPTS
            if kp is not None:
                (nkx, nky, nkz), _shift = kp
                kpts_use = (nkx, nky, nkz)

            input_data = {
                'control': {
                    'calculation': 'scf',
                    'prefix': 'ase_calc',
                    'pseudo_dir': os.path.abspath('.'),
                    'outdir': OUTDIR,
                    'tprnfor': True,
                    'tstress': False
                },
                'system': system_dict,
                'electrons': electrons_dict
            }
            # --- ここから 追加（Espresso を作る “直前”）---
            calc_dir = Path("vib_ads_vac")
            calc_dir.mkdir(parents=True, exist_ok=True)

            if use_environ:
                # ①で見つけた environ_src を calc.directory にも置く（上書きOK）
                shutil.copy2(environ_src, calc_dir / "environ.in")
            # --- ここまで 追加 ---

            self.atoms.calc = Espresso(
                pseudopotentials=PSEUDOS,
                input_data=input_data,
                kpts=kpts_use,          # ← ここを変更
                profile=profile,
                directory="vib_ads_vac"
            )
            vib = Vibrations(self.atoms, indices=mobile, delta=DELTA, name=VIB_NAME)
            if use_environ:
                candidates = [
                    Path("vib_ads_vac/environ.in"),
                    Path(f"{VIB_NAME}.eq/environ.in"),
                ]
                if not any(p.exists() for p in candidates):
                    raise FileNotFoundError(
                        "environ.in の事前配布に失敗しています。"
                        "vib_ads_vac/ または vib_ads.eq/ に environ.in がありません。"
                    )
            _ = sanitize_vib_cache(VIB_NAME, natoms=len(self.atoms))
            vib.run()
            freqs_cm = vib.get_frequencies()
            valid = [f for f in freqs_cm if f > LOW_CUT_CM]
            hc_eV_cm = 1.239841984e-4
            zpe_eV = 0.5 * hc_eV_cm * float(np.sum(valid))
            #S_vib(T)の算出
            T = 298.15  # ユーザー入力にしてもOK
            c_light = 2.99792458e10
            freqs_hz = np.array(valid) * c_light
            S_vib = 0.0
            for nu in freqs_hz:
                x = h * nu / (kB * T)
                S_vib += kB * (x / (np.exp(x) - 1.0) - np.log(1.0 - np.exp(-x)))
            S_vib_JmolK = S_vib * 96485.33212           
            vib.summary()
            end_time = datetime.now()
            # 経過時間を算出（分・秒）
            elapsed = end_time - start_time
            total_seconds = elapsed.total_seconds()
            minutes = int(total_seconds // 60)
            seconds = int(total_seconds % 60)
            ts = end_time.strftime('%Y%m%d-%H%M%S')
            base = f"{VIB_NAME}_{ts}"
            txt_path = Path(f"zpe_summary_{base}.txt")
            csv_path = Path(f"zpe_freqs_{base}.csv")
            with open(txt_path, 'w', encoding='utf-8') as f:
                f.write("# ZPE summary (ASE Vibrations)\n")
                f.write(f"calc_start_time: {start_time.strftime('%Y-%m-%d %H:%M:%S')}\n")
                f.write(f"calc_end_time: {end_time.strftime('%Y-%m-%d %H:%M:%S')}\n")
                f.write(f"elapsed: {minutes}m {seconds}s ({total_seconds:.1f} s)\n")
                f.write(f"qe_input: {self.path_var.get()}\n")
                f.write(f"selected_indices (mobile): {mobile}\n")
                f.write(f"delta (Å): {DELTA}\n")
                f.write(f"kpts: {kpts_use}\n")
                ecutwfc = system_dict.get('ecutwfc', 'N/A')
                ecutrho = system_dict.get('ecutrho', 'N/A')
                f.write(f"ecutwfc/ecutrho: {ecutwfc}/{ecutrho}\n")
                f.write(f"low_cut (cm^-1): {LOW_CUT_CM}\n")
                f.write(f"ZPE (eV): {zpe_eV:.6f}\n")
                f.write(f"S_vib(298.15 K): {S_vib_JmolK:.2f} J/mol/K\n")
                f.write(f"use_environ: {self.use_environ_var.get()}\n")
                f.write(f"calculation_mode: {'new' if new_calc else 'continue'}\n")
                f.write(f"cache_checked: {cleaned['checked_files']}\n")
                f.write(f"cache_deleted: {cleaned['deleted_files']}\n")
            with open(csv_path, 'w', newline='', encoding='utf-8') as f:
                writer = csv.writer(f)
                writer.writerow(["frequency_cm^-1"])
                for v in freqs_cm:
                    writer.writerow([f"{v:.6f}"])
            self.after(0, messagebox.showinfo, "完了",
                    (f"周波数 (cm^-1):\n{np.array2string(freqs_cm, precision=1, separator=', ')}\n\n"
                    f"ZPE (cut>{LOW_CUT_CM:.0f} cm^-1) = {zpe_eV:.6f} eV\n\n"
                    f"保存:\n  {txt_path}\n  {csv_path}\n\n"
                    f"経過時間: {minutes}m {seconds}s ({total_seconds:.1f} s)"))

        except Exception:
            err_text = traceback.format_exc()
            ts = datetime.now().strftime('%Y%m%d-%H%M%S')
            log_path = Path(f"vib_error_{ts}.log")
            log_path.write_text(err_text, encoding='utf-8')
            self.after(0, messagebox.showerror, "計算失敗", f"ログ: {log_path}\n\n{err_text}")
        finally:
            # ボタンを再び有効化（この版は btn_calc / btn_newcalc）
            self.after(0, self.btn_calc.config, {"state": tk.NORMAL})
            self.after(0, self.btn_newcalc.config, {"state": tk.NORMAL})
            self.after(0, self.status_var.set, "待機中。必要なら再度[計算]または[新規計算]を押してください。")

if __name__ == "__main__":
    app = ZPEGUI()
    app.mainloop()

