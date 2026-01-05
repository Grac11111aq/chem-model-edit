from __future__ import annotations
import sys, re
from dataclasses import dataclass
from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import List, Tuple, Optional 
from copy import deepcopy
import math 

from PyQt5.QtCore import Qt, QUrl
from PyQt5.QtWidgets import QComboBox, QLabel, QLineEdit, QFileDialog 
from PyQt5.QtGui import QGuiApplication
from PyQt5.QtWidgets import (
    QApplication, QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton,
    QTableWidget, QTableWidgetItem, QCheckBox, QInputDialog, QMessageBox,
)
from PyQt5.QtWebEngineWidgets import QWebEngineView
from PyQt5.QtWidgets import QSlider

# ---------------------------------------------------------------------------
# データクラス Atom: 原子一つ分の情報を保持する
# ---------------------------------------------------------------------------
@dataclass
class Atom:
    symbol: str
    x: float
    y: float
    z: float
    f1: str
    f2: str
    f3: str

    def pdb_line(self, idx: int, res_name="UNK", chain_id="A", res_seq=1) -> str:
        atom_name_field = f" {self.symbol.strip():<3}"
        if len(atom_name_field) > 4 : atom_name_field = atom_name_field[:4]
        element_field = f"{self.symbol.strip():>2}"
        return (
            f"ATOM  {idx:5d} {atom_name_field:4s}"
            f" {res_name:<3s}"
            f" {chain_id:1s}{res_seq:4d}"
            f"    {self.x:8.3f}{self.y:8.3f}{self.z:8.3f}"
            f"{1.00:6.2f}{0.00:6.2f}      "
            f"    {element_field:2s}"
        )

    def xyz_line(self) -> str:
        return f"{self.symbol} {self.x:.10f} {self.y:.10f} {self.z:.10f} {self.f1} {self.f2} {self.f3}".rstrip()

# ---------------------------------------------------------------------------
# XYZ形式のテキストブロックをパースするための正規表現と関数
# ---------------------------------------------------------------------------
XYZ_PATTERN = re.compile(
    r"^\s*(?P<sym>[A-Za-z]{1,2})\s+(?P<x>-?\d+\.\d+)\s+(?P<y>-?\d+\.\d+)\s+(?P<z>-?\d+\.\d+)"
    r"(?:\s+(?P<f1>\S+)\s+(?P<f2>\S+)\s+(?P<f3>\S+))?"
)

def parse_xyz_block(text: str) -> List[Atom]:
    atoms: List[Atom] = []
    for ln in text.strip().splitlines():
        m = XYZ_PATTERN.match(ln)
        if not m: continue
        atoms.append(
            Atom(
                m.group("sym"), float(m.group("x")), float(m.group("y")), float(m.group("z")),
                m.group("f1") or "", m.group("f2") or "", m.group("f3") or ""
            )
        )
    return atoms

# ---------------------------------------------------------------------------
# AtomTableクラス
# ---------------------------------------------------------------------------
class AtomTable(QTableWidget):
    COLS = ["Atom", "x", "y", "z", "flag1", "flag2", "flag3"]
    def __init__(self, name: str, main_window_ref: 'MainWindow'):
        super().__init__(0, len(self.COLS))
        self.setHorizontalHeaderLabels(self.COLS)
        self.setSelectionBehavior(self.SelectRows)
        self.setEditTriggers(self.DoubleClicked | self.SelectedClicked)
        self.name = name
        self.main_window_ref = main_window_ref
        self._original_atoms_on_load: List[Atom] = []
        self._currently_displayed_atoms: List[Atom] = []

    def _set_displayed_atoms(self, atoms_to_display: List[Atom]):
        self._currently_displayed_atoms = atoms_to_display
        self.setRowCount(0)
        for a in atoms_to_display: self._add_row(a)
        self.resizeColumnsToContents()

    def load_atoms(self, atoms: List[Atom]):
        self._original_atoms_on_load = deepcopy(atoms)
        self._set_displayed_atoms(deepcopy(atoms))

    def _add_row(self, a: Atom):
        r = self.rowCount(); self.insertRow(r)
        for c, txt in enumerate([a.symbol, f"{a.x:.10f}", f"{a.y:.10f}", f"{a.z:.10f}", a.f1, a.f2, a.f3]):
            it = QTableWidgetItem(txt)
            it.setTextAlignment(Qt.AlignRight|Qt.AlignVCenter)
            self.setItem(r,c,it)

    def to_atoms(self) -> List[Atom]:
        atoms_list = []
        for r in range(self.rowCount()):
            try:
                symbol = self.item(r,0).text(); x=float(self.item(r,1).text()); y=float(self.item(r,2).text()); z=float(self.item(r,3).text())
                f1=self.item(r,4).text() if self.item(r,4) else ""; f2=self.item(r,5).text() if self.item(r,5) else ""; f3=self.item(r,6).text() if self.item(r,6) else ""
                atoms_list.append(Atom(symbol,x,y,z,f1,f2,f3))
            except (AttributeError, ValueError): continue
        return atoms_list

    def copy_clip(self):
        lines = [
            f"{a.symbol:<2} {a.x:>16.10f} {a.y:>16.10f} {a.z:>16.10f} {a.f1} {a.f2} {a.f3}".rstrip()
            for a in self.to_atoms()
        ]
        QGuiApplication.clipboard().setText("\n".join(lines))

    def paste_clip(self):
        atoms = parse_xyz_block(QGuiApplication.clipboard().text())
        if atoms:
            self.load_atoms(atoms)
            self.main_window_ref._reset_shift_ui_for_table(self)
        else: QMessageBox.warning(self, "Paste", "クリップボードのXYZデータが無効です。")
        self.main_window_ref.refresh()

    def shift_selected(self, dx: float, dy: float, dz: float):
        coord_delta = {1:dx, 2:dy, 3:dz}
        for idx in self.selectedIndexes():
            col = idx.column()
            if col in coord_delta:
                try:
                    cur = float(self.item(idx.row(), col).text())
                    self.item(idx.row(), col).setText(f"{cur + coord_delta[col]:.10f}")
                except (ValueError, AttributeError): pass

    def sel_rows(self) -> List[int]: return sorted({i.row() for i in self.selectedIndexes()})

    def apply_global_shift(self, shift_x: float, shift_y: float, shift_z: float):
        if not self._original_atoms_on_load: return
        shifted_atoms = [Atom(a.symbol, a.x-shift_x, a.y-shift_y, a.z-shift_z, a.f1,a.f2,a.f3) for a in self._original_atoms_on_load]
        self._set_displayed_atoms(shifted_atoms)

    def revert_global_shift(self):
        if not self._original_atoms_on_load: self._set_displayed_atoms([]); return
        self._set_displayed_atoms(deepcopy(self._original_atoms_on_load))
    
    def commit_current_state(self):
        self._original_atoms_on_load = deepcopy(self.to_atoms())

# ---------------------------------------------------------------------------
# DistanceTableクラス
# ---------------------------------------------------------------------------
class DistanceTable(QTableWidget):
    COLS = ["Row #", "Atom A", "Atom B", "Distance (Å)"]
    def __init__(self):
        super().__init__(0, len(self.COLS))
        self.setHorizontalHeaderLabels(self.COLS); self.setEditTriggers(self.NoEditTriggers)
        self.setSelectionBehavior(self.SelectRows)

    def update_distances(self, atoms_a: List[Atom], atoms_b: List[Atom]):
        self.setRowCount(0); min_rows = min(len(atoms_a), len(atoms_b))
        for i in range(min_rows):
            a1,a2 = atoms_a[i], atoms_b[i]
            try: distance = math.sqrt((a1.x-a2.x)**2 + (a1.y-a2.y)**2 + (a1.z-a2.z)**2)
            except AttributeError: distance = float('nan')
            r = self.rowCount(); self.insertRow(r)
            for c,val in enumerate([str(i+1), a1.symbol, a2.symbol, f"{distance:.4f}" if not math.isnan(distance) else "Error"]):
                item = QTableWidgetItem(val); item.setTextAlignment(Qt.AlignRight|Qt.AlignVCenter if c in [0,3] else Qt.AlignVCenter)
                self.setItem(r,c,item)
        self.resizeColumnsToContents()

    def get_formatted_text(self) -> str:
        if self.rowCount()==0: return ""
        lines = ["\t".join(self.COLS)] + [
            "\t".join(self.item(r,c).text() if self.item(r,c) else "" for c in range(self.columnCount()))
            for r in range(self.rowCount())
        ]
        return "\n".join(lines)

# ---------------------------------------------------------------------------
# HTMLとJavaScript (NGL Viewer用)
# ---------------------------------------------------------------------------
HTML_TEMPLATE_FOR_APP="""<!DOCTYPE html><html><head><meta charset='utf-8'><script src='https://unpkg.com/ngl@2/dist/ngl.js'></script></head>
<body style='margin:0'><div id='viewport' style='width:100%;height:100vh'></div>
<script>
var stage=new NGL.Stage('viewport');
function clearByName(name){stage.compList.slice().forEach(c=>{if(c.name===name)stage.removeComponent(c);});}
function loadPDB(name,str,colorVal,op,showLabel){
 clearByName(name);
 var blob=new Blob([str],{type:'text/plain'});
  stage.loadFile(blob,{ext:'pdb',name:name}).then(function(c){
    var repParams = {opacity:op};
    if (colorVal && colorVal.trim() !== "") { // colorValがnullや空文字列でないことを確認
        repParams.color = colorVal;
    }
    c.addRepresentation('ball+stick', repParams);
    if(showLabel){
        c.addRepresentation('label',{labelType:'serial', color:'blue', scale:0.8});
    }
    stage.autoView(); 
 });
}
</script></body></html>"""

# ---------------------------------------------------------------------------
# Viewerクラス
# ---------------------------------------------------------------------------
class Viewer(QWebEngineView):
    def __init__(self):
        super().__init__()
        tmp_html_file = NamedTemporaryFile(delete=False, suffix='.html', mode='w', encoding='utf-8')
        tmp_html_file.write(HTML_TEMPLATE_FOR_APP) 
        tmp_html_file.close()
        self.load(QUrl.fromLocalFile(tmp_html_file.name))
        self._ready = False
        self.loadFinished.connect(lambda _: setattr(self, '_ready', True))

    def update(self,name:str,pdb:str,color_val:Optional[str],op:float,show_label:bool):
        if self._ready:
            js_color_arg = f"'{color_val}'" if color_val is not None else "null"
            self.page().runJavaScript(
                f"loadPDB('{name}', `{pdb}`, {js_color_arg}, {op}, {str(show_label).lower()});")

    def remove(self,name:str):
        if self._ready: self.page().runJavaScript(f"clearByName('{name}');")

    def set_background(self, color: str):
        if self._ready: self.page().runJavaScript(f"stage.setParameters({{ backgroundColor: '{color}' }});")

# ---------------------------------------------------------------------------
# MainWindowクラス
# ---------------------------------------------------------------------------
class MainWindow(QWidget):
    def __init__(self):
        super().__init__(); self.setWindowTitle('Atomic Coordinate Editor'); self.resize(1600,820)
        
        self.tA = AtomTable('A', self); self.tB = AtomTable('B', self)
        self.tC = DistanceTable(); self.view = Viewer()
        
        # --- グローバルコントロールUI要素 ---
        self.bg_label = QLabel("背景色:")
        self.bg_combo = QComboBox()
        self.bg_combo.addItems(["White", "Black"])
        self.bg_combo.currentIndexChanged.connect(self.refresh)
        
        # --- Table A/B 共通表示設定 (これらは refresh で直接値が読まれる) ---
        self.chkA = QCheckBox('表示(A)')
        self.trA = QCheckBox('半透明(A)') # この半透明チェックボックスは現在 refresh メソッドで直接は使用されていない
        self.chkB = QCheckBox('表示(B)')
        self.trB = QCheckBox('半透明(B)') # この半透明チェックボックスは現在 refresh メソッドで直接は使用されていない
        for c in (self.chkA, self.trA, self.chkB, self.trB): 
            c.stateChanged.connect(self.refresh) # 半透明のチェックボックスはrefreshに影響を与えるが、個別の処理はない

        # --- 各テーブル固有のコントロール (スライダー、ラベル表示、色モードなど) ---
        # Table A
        self.opacityA_label = QLabel("不透明度 A:") # ****** 属性をここで定義 ******
        self.sliderA = QSlider(Qt.Horizontal); self.sliderA.setRange(0,100); self.sliderA.setValue(100)
        self.sliderA.valueChanged.connect(self.refresh)
        self.chkLabelA = QCheckBox("番号表示(A)")
        self.chkLabelA.stateChanged.connect(self.refresh)
        self.colorModeA_label = QLabel("色モード (A):")
        self.colorModeA_combo = QComboBox()
        self.colorModeA_combo.addItems(["固定色 (Skyblue)", "デフォルト (CPK)", "RGB/名称指定"])
        self.colorModeA_combo.currentIndexChanged.connect(lambda: self._on_color_mode_changed(self.colorModeA_combo, self.rgbInputA_label, self.rgbInputA_edit))
        self.rgbInputA_label = QLabel("RGB/名称 (A):")
        self.rgbInputA_edit = QLineEdit("skyblue") 
        self.rgbInputA_label.setVisible(False)
        self.rgbInputA_edit.setVisible(False)
        self.rgbInputA_edit.editingFinished.connect(self.refresh)

        # Table B
        self.opacityB_label = QLabel("不透明度 B:") # ****** 属性をここで定義 ******
        self.sliderB = QSlider(Qt.Horizontal); self.sliderB.setRange(0,100); self.sliderB.setValue(100)
        self.sliderB.valueChanged.connect(self.refresh)
        self.chkLabelB = QCheckBox("番号表示(B)")
        self.chkLabelB.stateChanged.connect(self.refresh)
        self.colorModeB_label = QLabel("色モード (B):")
        self.colorModeB_combo = QComboBox()
        self.colorModeB_combo.addItems(["固定色 (Salmon)", "デフォルト (CPK)", "RGB/名称指定"])
        self.colorModeB_combo.currentIndexChanged.connect(lambda: self._on_color_mode_changed(self.colorModeB_combo, self.rgbInputB_label, self.rgbInputB_edit))
        self.rgbInputB_label = QLabel("RGB/名称 (B):")
        self.rgbInputB_edit = QLineEdit("salmon") 
        self.rgbInputB_label.setVisible(False)
        self.rgbInputB_edit.setVisible(False)
        self.rgbInputB_edit.editingFinished.connect(self.refresh)

        # --- ボタン類 (Table A, B) ---
        pA,cA,sA=QPushButton('Paste A'),QPushButton('Copy A'),QPushButton('Shift Selected Cells A')
        pB,cB,sB=QPushButton('Paste B'),QPushButton('Copy B'),QPushButton('Shift Selected Cells B')
        align=QPushButton('Align → B')
        commitA = QPushButton('現在の座標を確定 (A)'); commitB = QPushButton('現在の座標を確定 (B)')
        
        pA.clicked.connect(self.tA.paste_clip); pB.clicked.connect(self.tB.paste_clip)
        cA.clicked.connect(self.tA.copy_clip); cB.clicked.connect(self.tB.copy_clip)
        sA.clicked.connect(lambda:self.shift_dialog(self.tA)); sB.clicked.connect(lambda:self.shift_dialog(self.tB))
        align.clicked.connect(self.align_copy)
        commitA.clicked.connect(lambda: self._commit_table_state(self.tA))
        commitB.clicked.connect(lambda: self._commit_table_state(self.tB))

        # --- 選択原子座標表示 & グローバルシフトUI ---
        self.coordA_x_display = QLineEdit(); self.coordA_x_display.setReadOnly(True); self.coordA_x_display.setPlaceholderText("X")
        self.coordA_y_display = QLineEdit(); self.coordA_y_display.setReadOnly(True); self.coordA_y_display.setPlaceholderText("Y")
        self.coordA_z_display = QLineEdit(); self.coordA_z_display.setReadOnly(True); self.coordA_z_display.setPlaceholderText("Z")
        self.chkShiftTableA = QCheckBox("Shift Table A by Selected Coords")

        self.coordB_x_display = QLineEdit(); self.coordB_x_display.setReadOnly(True); self.coordB_x_display.setPlaceholderText("X")
        self.coordB_y_display = QLineEdit(); self.coordB_y_display.setReadOnly(True); self.coordB_y_display.setPlaceholderText("Y")
        self.coordB_z_display = QLineEdit(); self.coordB_z_display.setReadOnly(True); self.coordB_z_display.setPlaceholderText("Z")
        self.chkShiftTableB = QCheckBox("Shift Table B by Selected Coords")

        self.tA.itemSelectionChanged.connect(self._update_coord_display_A)
        self.tB.itemSelectionChanged.connect(self._update_coord_display_B)
        self.chkShiftTableA.stateChanged.connect(self._toggle_table_shift_A)
        self.chkShiftTableB.stateChanged.connect(self._toggle_table_shift_B)

        # --- レイアウト設定 ---
        main_layout = QVBoxLayout(self) # メインの垂直レイアウト

        # 上部全体コントロール (背景色のみ)
        global_controls_layout = QHBoxLayout()
        global_controls_layout.addWidget(self.bg_label); global_controls_layout.addWidget(self.bg_combo)
        global_controls_layout.addStretch(1) # 右に寄せるためのスペーサー
        main_layout.addLayout(global_controls_layout)
        
        # テーブルごとの詳細コントロールセクション (水平レイアウト)
        table_specific_controls_layout = QHBoxLayout()
        
        # Table A コントロール群
        #colA_specific_controls = QVBoxLayout()

        #table_specific_controls_layout.addLayout(colA_specific_controls)

        # Table B コントロール群
        #colB_specific_controls = QVBoxLayout()
        #table_specific_controls_layout.addLayout(colB_specific_controls)
        
        #main_layout.addLayout(table_specific_controls_layout) # メインレイアウトに追加
        
        # テーブル群のレイアウト
        colA_layout, colB_layout, colC_layout = QVBoxLayout(), QVBoxLayout(), QVBoxLayout()
        
        # Table Aカラム
        colA_layout.addWidget(QLabel(self.tA.name)); colA_layout.addWidget(self.tA,1)
        btn_layout_A = QHBoxLayout(); btn_layout_A.addWidget(pA); btn_layout_A.addWidget(cA); btn_layout_A.addWidget(sA)
        colA_layout.addLayout(btn_layout_A); colA_layout.addWidget(commitA)
        coordA_label = QLabel("Selected Atom A:"); colA_layout.addWidget(coordA_label)
        coordA_hlayout = QHBoxLayout()
        coordA_hlayout.addWidget(QLabel("X:")); coordA_hlayout.addWidget(self.coordA_x_display)
        coordA_hlayout.addWidget(QLabel("Y:")); coordA_hlayout.addWidget(self.coordA_y_display)
        coordA_hlayout.addWidget(QLabel("Z:")); coordA_hlayout.addWidget(self.coordA_z_display)
        colA_layout.addLayout(coordA_hlayout); colA_layout.addWidget(self.chkShiftTableA)
        # chkA と trA は各テーブルの列の下部に配置
        display_chk_layout_A = QVBoxLayout(); display_chk_layout_A.addWidget(self.chkA); display_chk_layout_A.addWidget(self.trA)
        display_chk_layout_A.addWidget(self.opacityA_label); display_chk_layout_A.addWidget(self.sliderA)
        display_chk_layout_A.addWidget(self.chkLabelA)
        display_chk_layout_A.addWidget(self.colorModeA_label); display_chk_layout_A.addWidget(self.colorModeA_combo)
        display_chk_layout_A.addWidget(self.rgbInputA_label); display_chk_layout_A.addWidget(self.rgbInputA_edit)
        display_chk_layout_A.addStretch(1) # 下に寄せる
        colA_layout.addLayout(display_chk_layout_A)

        # Table Bカラム
        colB_layout.addWidget(QLabel(self.tB.name)); colB_layout.addWidget(self.tB,1)
        btn_layout_B = QHBoxLayout(); btn_layout_B.addWidget(pB); btn_layout_B.addWidget(cB); btn_layout_B.addWidget(sB)
        colB_layout.addLayout(btn_layout_B); colB_layout.addWidget(commitB)
        coordB_label = QLabel("Selected Atom B:"); colB_layout.addWidget(coordB_label)
        coordB_hlayout = QHBoxLayout()
        coordB_hlayout.addWidget(QLabel("X:")); coordB_hlayout.addWidget(self.coordB_x_display)
        coordB_hlayout.addWidget(QLabel("Y:")); coordB_hlayout.addWidget(self.coordB_y_display)
        coordB_hlayout.addWidget(QLabel("Z:")); coordB_hlayout.addWidget(self.coordB_z_display)
        colB_layout.addLayout(coordB_hlayout); colB_layout.addWidget(self.chkShiftTableB)
        # chkB と trB は各テーブルの列の下部に配置
        display_chk_layout_B = QVBoxLayout(); display_chk_layout_B.addWidget(self.chkB); display_chk_layout_B.addWidget(self.trB)
        display_chk_layout_B.addWidget(self.opacityB_label); display_chk_layout_B.addWidget(self.sliderB)
        display_chk_layout_B.addWidget(self.chkLabelB)
        display_chk_layout_B.addWidget(self.colorModeB_label); display_chk_layout_B.addWidget(self.colorModeB_combo)
        display_chk_layout_B.addWidget(self.rgbInputB_label); display_chk_layout_B.addWidget(self.rgbInputB_edit)
        display_chk_layout_B.addStretch(1) # 下に寄せる
        colB_layout.addLayout(display_chk_layout_B); colB_layout.addWidget(align)

        # Table Cカラム
        colC_layout.addWidget(QLabel("Table C: Distances A-B")); colC_layout.addWidget(self.tC, 1)
        self.btnCopyC = QPushButton("クリップボードにコピー (Table C)")
        colC_layout.addWidget(self.btnCopyC); self.btnCopyC.clicked.connect(self._copy_table_c_to_clipboard)

        tables_layout=QHBoxLayout()
        tables_layout.addLayout(colA_layout,1); tables_layout.addLayout(colB_layout,1); tables_layout.addLayout(colC_layout,1)
        
        export_btn = QPushButton("HTML形式でエクスポート"); export_btn.clicked.connect(self._export_to_html)
        
        root_layout=QHBoxLayout(); root_layout.addLayout(tables_layout,3); root_layout.addWidget(self.view,4) # Viewerの比率を調整
        
        main_layout.addLayout(root_layout, 1); main_layout.addWidget(export_btn)
     
    def _on_color_mode_changed(self, combo_box: QComboBox, rgb_label: QLabel, rgb_edit: QLineEdit):
        is_rgb_mode = ("RGB" in combo_box.currentText())
        rgb_label.setVisible(is_rgb_mode)
        rgb_edit.setVisible(is_rgb_mode)
        self.refresh()

    def _get_color_for_table(self, table_char: str) -> Optional[str]:
        if table_char == 'A':
            mode_combo, fixed_color, rgb_edit, default_fixed_color = \
                self.colorModeA_combo, "skyblue", self.rgbInputA_edit, "skyblue"
        else:
            mode_combo, fixed_color, rgb_edit, default_fixed_color = \
                self.colorModeB_combo, "salmon", self.rgbInputB_edit, "salmon"

        selected_mode = mode_combo.currentText()

        if "デフォルト (CPK)" in selected_mode: return None
        elif "固定色" in selected_mode: return fixed_color
        elif "RGB/名称指定" in selected_mode:
            rgb_text = rgb_edit.text().strip()
            if rgb_text and (re.match(r"^#(?:[0-9a-fA-F]{3}){1,2}$", rgb_text) or \
                             re.match(r"^#(?:[0-9a-fA-F]{4}){1,2}$", rgb_text) or \
                             re.match(r"^[a-zA-Z]+$", rgb_text)):
                return rgb_text
            else:
                # QMessageBox.warning(self, "色指定エラー", f"Table {table_char} のRGB/名称指定「{rgb_text}」が無効です。固定色 ({default_fixed_color}) を使用します。")
                # rgb_edit.setText(default_fixed_color) # ユーザー入力を勝手に変更しない方が良いかもしれない
                return default_fixed_color # 無効な場合は固定色をフォールバック
        return None

    def refresh(self):
        bg = self.bg_combo.currentText().lower(); self.view.set_background(bg)
        atoms_A_current = self.tA.to_atoms(); atoms_B_current = self.tB.to_atoms()

        color_a = self._get_color_for_table('A')
        color_b = self._get_color_for_table('B')
        
        # trA, trB (半透明チェックボックス) の状態を opacity に反映する
        opacity_a = 0.3 if self.trA.isChecked() else self.sliderA.value()/100.0
        opacity_b = 0.3 if self.trB.isChecked() else self.sliderB.value()/100.0


        if self.chkA.isChecked() and atoms_A_current:
            pdbA = self._atoms_to_pdb(atoms_A_current)
            self.view.update('A', pdbA, color_a, opacity_a, self.chkLabelA.isChecked())
        else: self.view.remove('A')

        if self.chkB.isChecked() and atoms_B_current:
            pdbB = self._atoms_to_pdb(atoms_B_current)
            self.view.update('B', pdbB, color_b, opacity_b, self.chkLabelB.isChecked())
        else: self.view.remove('B')
        self._update_distance_table()
    
    @staticmethod
    def _atoms_to_pdb(atoms: List[Atom]) -> str:
        lines=["HEADER    GENERATED BY ATOM_GUI"] + [a.pdb_line(i+1) for i,a in enumerate(atoms)] + ["END"]
        return "\n".join(lines)
    
    @staticmethod
    def _atoms_to_pdb_for_file(atoms: List[Atom]) -> str:
        lines=["HEADER    GENERATED BY ATOM_GUI"] + [a.pdb_line(i+1) for i,a in enumerate(atoms)] + ["END"]
        return "\n".join(lines)

    def shift_dialog(self,tbl:AtomTable):
        txt,ok=QInputDialog.getText(self,'Shift Selected Cells','dx dy dz (Å) で入力 (スペースまたはカンマ区切り):');
        if not ok or not txt: return
        try: dx,dy,dz=map(float,re.split(r"[ ,]+",txt.strip()))
        except ValueError: QMessageBox.warning(self,'Error','3つの数値を正しく入力してください。'); return
        tbl.shift_selected(dx,dy,dz); self.refresh()

    def align_copy(self):
        selA_rows, selB_rows = self.tA.sel_rows(), self.tB.sel_rows()
        if not selA_rows: QMessageBox.warning(self,'Align','Table A で行を選択してください.'); return
        if len(selB_rows) != 1: QMessageBox.warning(self,'Align','Table B で基準となる行を1行だけ選択してください.'); return
        refB_idx = selB_rows[0]
        atomsA = self.tA.to_atoms(); atomsB = self.tB.to_atoms()
        if not (0 <= refB_idx < len(atomsB)): QMessageBox.warning(self,'Align',f'Table B で選択された基準行 {refB_idx+1} は無効です。'); return
        refA_idx_for_delta = refB_idx 
        if not (0 <= refA_idx_for_delta < len(atomsA)): QMessageBox.warning(self,'Align',f'Table A の対応する基準行 {refA_idx_for_delta+1} は無効です。'); return
        ref_atom_A = atomsA[refA_idx_for_delta]; ref_atom_B = atomsB[refB_idx]
        dx = ref_atom_B.x-ref_atom_A.x; dy = ref_atom_B.y-ref_atom_A.y; dz = ref_atom_B.z-ref_atom_A.z
        for rA_idx in selA_rows:
            if rA_idx < len(atomsA): 
                atom_from_A = atomsA[rA_idx]; rB_target_idx = rA_idx 
                if rB_target_idx < self.tB.rowCount():
                    data_to_set = [atom_from_A.symbol, f"{atom_from_A.x+dx:.10f}", f"{atom_from_A.y+dy:.10f}", f"{atom_from_A.z+dz:.10f}", atom_from_A.f1, atom_from_A.f2, atom_from_A.f3]
                    for c_idx, val_str in enumerate(data_to_set):
                        item = self.tB.item(rB_target_idx, c_idx)
                        if item is None: item = QTableWidgetItem(); self.tB.setItem(rB_target_idx, c_idx, item)
                        item.setText(val_str)
                        if 0 < c_idx < 4: item.setTextAlignment(Qt.AlignRight|Qt.AlignVCenter)
        self.refresh()

    def _update_coord_display(self, table: AtomTable, x_edit: QLineEdit, y_edit: QLineEdit, z_edit: QLineEdit):
        sel_rows = table.sel_rows()
        if sel_rows and 0 <= sel_rows[0] < table.rowCount():
            row = sel_rows[0]
            try: x_edit.setText(table.item(row,1).text()); y_edit.setText(table.item(row,2).text()); z_edit.setText(table.item(row,3).text())
            except AttributeError: x_edit.clear(); y_edit.clear(); z_edit.clear()
        else: x_edit.clear(); y_edit.clear(); z_edit.clear()

    def _update_coord_display_A(self): self._update_coord_display(self.tA, self.coordA_x_display, self.coordA_y_display, self.coordA_z_display)
    def _update_coord_display_B(self): self._update_coord_display(self.tB, self.coordB_x_display, self.coordB_y_display, self.coordB_z_display)

    def _toggle_table_shift(self, state: int, table: AtomTable, x_edit: QLineEdit, y_edit: QLineEdit, z_edit: QLineEdit, chkBox: QCheckBox):
        if state == Qt.Checked:
            if not table._original_atoms_on_load: QMessageBox.warning(self,"シフトエラー",f"{table.name}に元のデータ無"); chkBox.setChecked(False); return
            try: shift_x=float(x_edit.text()); shift_y=float(y_edit.text()); shift_z=float(z_edit.text())
            except ValueError: QMessageBox.warning(self,"シフトエラー","X,Y,Z座標が無効"); chkBox.setChecked(False); return
            table.apply_global_shift(shift_x, shift_y, shift_z)
        else: table.revert_global_shift()
        if table == self.tA: self._update_coord_display_A()
        else: self._update_coord_display_B()
        self.refresh()

    def _toggle_table_shift_A(self, state): self._toggle_table_shift(state, self.tA, self.coordA_x_display, self.coordA_y_display, self.coordA_z_display, self.chkShiftTableA)
    def _toggle_table_shift_B(self, state): self._toggle_table_shift(state, self.tB, self.coordB_x_display, self.coordB_y_display, self.coordB_z_display, self.chkShiftTableB)

    def _reset_shift_ui_for_table(self, table: AtomTable):
        chk,xd,yd,zd = (self.chkShiftTableA,self.coordA_x_display,self.coordA_y_display,self.coordA_z_display) if table==self.tA else (self.chkShiftTableB,self.coordB_x_display,self.coordB_y_display,self.coordB_z_display)
        if chk.isChecked(): chk.blockSignals(True); chk.setChecked(False); chk.blockSignals(False)
        xd.clear(); yd.clear(); zd.clear()

    def _update_distance_table(self): self.tC.update_distances(self.tA.to_atoms(), self.tB.to_atoms())

    def _copy_table_c_to_clipboard(self):
        txt=self.tC.get_formatted_text()
        if not txt: QMessageBox.information(self,"コピー (Table C)","Table C は空です。"); return
        QGuiApplication.clipboard().setText(txt); QMessageBox.information(self,"コピー (Table C)","Table C の内容をコピーしました。")

    def _commit_table_state(self, table: AtomTable):
        if table.rowCount()==0: QMessageBox.warning(self,"座標の確定",f"Table {table.name} は空です。"); return
        reply = QMessageBox.question(self,'座標の確定',f"Table {table.name} の現在の表示座標を新しい基準として確定しますか？\n\nこの操作は元に戻せません。",QMessageBox.Yes|QMessageBox.No,QMessageBox.No)
        if reply==QMessageBox.Yes: table.commit_current_state(); self._reset_shift_ui_for_table(table); QMessageBox.information(self,"成功",f"Table {table.name} の座標が確定されました。")

    def _export_to_html(self):
        filePath, _ = QFileDialog.getSaveFileName(self, "HTML形式で保存", "", "HTML Files (*.html);;All Files (*)")
        if not filePath: return

        is_a_visible = self.chkA.isChecked() and self.tA.rowCount() > 0
        is_b_visible = self.chkB.isChecked() and self.tB.rowCount() > 0
        if not is_a_visible and not is_b_visible: QMessageBox.information(self, "エクスポート", "表示されている3Dモデルがありません。"); return
        
        bg_color = self.bg_combo.currentText().lower()
        
        js_data_definitions = []
        js_promise_creations = []
        js_promise_names = []

        if is_a_visible:
            pdb_a_str = self._atoms_to_pdb_for_file(self.tA.to_atoms())
            opacity_a = 0.3 if self.trA.isChecked() else self.sliderA.value() / 100.0 # 半透明チェックを反映
            label_a = str(self.chkLabelA.isChecked()).lower()
            color_a_val = self._get_color_for_table('A')
            js_color_a_param = f"'{color_a_val}'" if color_a_val else "null"
            
            js_data_definitions.append(f"const pdbDataA = `{pdb_a_str}`;")
            js_data_definitions.append(f"const blobA = new Blob([pdbDataA], {{ type: 'text/plain' }});")
            js_promise_creations.append(f"""const promiseA = stage.loadFile(blobA, {{ ext: 'pdb', name: 'ModelA' }}).then(function (c) {{
                var repParamsA = {{ opacity: {opacity_a} }};
                if ({js_color_a_param} !== null) {{ repParamsA.color = {js_color_a_param}; }}
                c.addRepresentation('ball+stick', repParamsA);
                if ({label_a}) {{
                    c.addRepresentation('label', {{ labelType: 'serial', color: 'blue', scale: 0.8 }});
                }}
            }});""")
            js_promise_names.append("promiseA")

        if is_b_visible:
            pdb_b_str = self._atoms_to_pdb_for_file(self.tB.to_atoms())
            opacity_b = 0.3 if self.trB.isChecked() else self.sliderB.value() / 100.0 # 半透明チェックを反映
            label_b = str(self.chkLabelB.isChecked()).lower()
            color_b_val = self._get_color_for_table('B')
            js_color_b_param = f"'{color_b_val}'" if color_b_val else "null"
            
            js_data_definitions.append(f"const pdbDataB = `{pdb_b_str}`;")
            js_data_definitions.append(f"const blobB = new Blob([pdbDataB], {{ type: 'text/plain' }});")
            js_promise_creations.append(f"""const promiseB = stage.loadFile(blobB, {{ ext: 'pdb', name: 'ModelB' }}).then(function (c) {{
                var repParamsB = {{ opacity: {opacity_b} }};
                if ({js_color_b_param} !== null) {{ repParamsB.color = {js_color_b_param}; }}
                c.addRepresentation('ball+stick', repParamsB);
                if ({label_b}) {{
                    c.addRepresentation('label', {{ labelType: 'serial', color: 'blue', scale: 0.8 }});
                }}
            }});""")
            js_promise_names.append("promiseB")
        
        ngl_cdn_url = "https://unpkg.com/ngl@2/dist/ngl.js"

        script_content = f"""
        document.addEventListener('DOMContentLoaded', function () {{
            try {{ 
                if (typeof NGL === 'undefined' || !NGL.Stage) {{
                    alert('NGLライブラリがロードされていません。インターネット接続を確認するか、HTML内のCDN URLを確認してください。');
                    console.error('NGL library not loaded or NGL.Stage is undefined.');
                    return;
                }}
                var stage = new NGL.Stage('viewport');
                stage.setParameters({{ backgroundColor: '{bg_color}' }});
                {'\n'.join(js_data_definitions)} 
                {'\n'.join(js_promise_creations)}
                Promise.all([{', '.join(js_promise_names)}]).then(function() {{
                    stage.autoView();
                }}).catch(function(error) {{
                    console.error("Error during Promise.all (loading models or autoView): ", error);
                    alert("モデルの読み込みまたは表示調整中にエラーが発生しました。詳細はブラウザのコンソールを確認してください。");
                }});
            }} catch (e) {{
                console.error("General JavaScript error during NGL setup: ", e);
                alert("NGLビューアの初期化中に予期せぬエラーが発生しました: " + e.message + "\\n詳細はブラウザのコンソールを確認してください。");
            }}
        }});
"""
        
        html_template = f"""<!DOCTYPE html><html lang="ja"><head><meta charset="utf-8"><title>NGL View Export</title></head>
<body style="margin:0; overflow:hidden;"><script src="{ngl_cdn_url}"></script><div id="viewport" style="width:100vw; height:100vh;"></div>
<script>{script_content}</script></body></html>"""

        try:
            with open(filePath, 'w', encoding='utf-8') as f: f.write(html_template)
            QMessageBox.information(self, "成功", f"HTMLファイルが正常に保存されました:\n{filePath}")
        except Exception as e: QMessageBox.critical(self, "エラー", f"ファイルの保存中にエラーが発生しました:\n{e}")

# ---------------------------------------------------------------------------
# アプリケーションの起動
# ---------------------------------------------------------------------------
if __name__=='__main__':
     app = QApplication(sys.argv); win = MainWindow(); win.show(); sys.exit(app.exec_())