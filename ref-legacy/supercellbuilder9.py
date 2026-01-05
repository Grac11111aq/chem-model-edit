# coding: utf-8
from __future__ import annotations
import sys, re
from dataclasses import dataclass
from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import List, Tuple, Optional, Dict
from copy import deepcopy
import math

from PyQt5.QtCore import Qt, QUrl
from PyQt5.QtWidgets import QGroupBox, QGridLayout, QFrame
from PyQt5.QtWidgets import QComboBox, QLabel, QLineEdit, QFileDialog
from PyQt5.QtGui import QGuiApplication
from PyQt5.QtWidgets import (
    QApplication, QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton,
    QTableWidget, QTableWidgetItem, QCheckBox, QInputDialog, QMessageBox,
)
from PyQt5.QtWebEngineWidgets import QWebEngineView
from PyQt5.QtWidgets import QSlider

# (Atom, parse_xyz_block, AtomTable, Viewer, HTML_TEMPLATE の各クラス/定義は変更ありません)
@dataclass
class Atom:
    symbol: str; x: float; y: float; z: float; f1: str; f2: str; f3: str
    def get_shifted(self, shift_vec: Tuple[float, float, float]) -> Atom:
        return Atom(self.symbol, self.x + shift_vec[0], self.y + shift_vec[1], self.z + shift_vec[2], self.f1, self.f2, self.f3)
    def pdb_line(self, idx: int, res_name="UNK", chain_id="A", res_seq=1) -> str:
        atom_name_field = f" {self.symbol.strip():<3}"
        if len(atom_name_field) > 4: atom_name_field = atom_name_field[:4]
        element_field = f"{self.symbol.strip():>2}"
        return (f"ATOM  {idx:5d} {atom_name_field:4s} {res_name:<3s} {chain_id:1s}{res_seq:4d}    "
                f"{self.x:8.3f}{self.y:8.3f}{self.z:8.3f}{1.00:6.2f}{0.00:6.2f}          {element_field:2s}")
    def xyz_line(self) -> str:
        return f"{self.symbol:<3s} {self.x:15.10f} {self.y:15.10f} {self.z:15.10f}  {self.f1} {self.f2} {self.f3}".rstrip()

def parse_xyz_block(text: str) -> List[Atom]:
    atoms: List[Atom] = []
    pattern = re.compile(
        r"^\s*(?P<sym>[A-Za-z]{1,3})\s+"
        r"(?P<x>-?\d+\.?\d*(?:[eE][+-]?\d+)?)\s+"
        r"(?P<y>-?\d+\.?\d*(?:[eE][+-]?\d+)?)\s+"
        r"(?P<z>-?\d+\.?\d*(?:[eE][+-]?\d+)?)\s*"
        r"(?P<f1>\S*)\s*(?P<f2>\S*)\s*(?P<f3>\S*)\s*$"
    )
    lines = [line for line in text.strip().splitlines() if not line.strip().startswith('#')]
    for line in lines:
        match = pattern.match(line)
        if match:
            data = match.groupdict()
            atoms.append(Atom(symbol=data['sym'], x=float(data['x']), y=float(data['y']), z=float(data['z']),
                              f1=data['f1'], f2=data['f2'], f3=data['f3']))
    return atoms

class AtomTable(QTableWidget):
    COLS = ["Atom", "x", "y", "z", "flag1", "flag2", "flag3"]
    def __init__(self, name: str, main_window_ref: 'MainWindow'):
        super().__init__(0, len(self.COLS)); self.setHorizontalHeaderLabels(self.COLS)
        self.setSelectionBehavior(self.SelectRows); self.setEditTriggers(self.DoubleClicked | self.SelectedClicked)
        self.name = name; self.main_window_ref = main_window_ref
    def load_atoms(self, atoms: List[Atom]):
        self.setRowCount(0)
        for a in atoms: self._add_row(a)
        self.resizeColumnsToContents()
    def _add_row(self, a: Atom):
        r = self.rowCount(); self.insertRow(r)
        row_data = [a.symbol, f"{a.x:.10f}", f"{a.y:.10f}", f"{a.z:.10f}", a.f1, a.f2, a.f3]
        for c, txt in enumerate(row_data):
            it = QTableWidgetItem(txt)
            if c > 0 and c < 4: it.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
            self.setItem(r,c,it)
    def to_atoms(self) -> List[Atom]:
        atoms_list = []
        for r in range(self.rowCount()):
            try:
                atoms_list.append(Atom(
                    symbol=self.item(r,0).text(), x=float(self.item(r,1).text()),
                    y=float(self.item(r,2).text()), z=float(self.item(r,3).text()),
                    f1=self.item(r,4).text() if self.item(r,4) and self.item(r,4).text() else "",
                    f2=self.item(r,5).text() if self.item(r,5) and self.item(r,5).text() else "",
                    f3=self.item(r,6).text() if self.item(r,6) and self.item(r,6).text() else ""
                ))
            except (AttributeError, ValueError): continue
        return atoms_list
    def copy_clip(self):
        atoms = self.to_atoms()
        if not atoms:
            QMessageBox.information(self, "コピー", f"Table {self.name} は空です。"); return
        if self.name == 'C (Supercell)':
            try:
                main_win = self.main_window_ref; na = main_win.supercell_na; nb = main_win.supercell_nb
                a_orig = float(main_win.lattice_edits['a'].text()); b_orig = float(main_win.lattice_edits['b'].text())
                c_orig = float(main_win.lattice_edits['c'].text()); cosab_orig = float(main_win.lattice_edits['cosab'].text())
                cosac_orig = float(main_win.lattice_edits['cosac'].text()); cosbc_orig = float(main_win.lattice_edits['cosbc'].text())
                va_orig, vb_orig, vc_orig = main_win.supercell_va_orig, main_win.supercell_vb_orig, main_win.supercell_vc_orig
                a_new, b_new, c_new = a_orig * na, b_orig * nb, c_orig
                va_new = tuple(c * na for c in va_orig); vb_new = tuple(c * nb for c in vb_orig); vc_new = vc_orig
                param_lines = [
                    "# Supercell Lattice Parameters (a, b, c, angles)", f"a     = {a_new:.8e}", f"b     = {b_new:.8e}",
                    f"c     = {c_new:.8e}", f"cosab = {cosab_orig:.8e}", f"cosac = {cosac_orig:.8e}", f"cosbc = {cosbc_orig:.8e}"
                ]
                qe_cell_lines = [
                    "CELL_PARAMETERS { angstrom }", f"   {va_new[0]:15.10f} {va_new[1]:15.10f} {va_new[2]:15.10f}",
                    f"   {vb_new[0]:15.10f} {vb_new[1]:15.10f} {vb_new[2]:15.10f}", f"   {vc_new[0]:15.10f} {vc_new[1]:15.10f} {vc_new[2]:15.10f}",
                ]
                atom_lines = ["# Atom Coordinates"] + [a.xyz_line() for a in atoms]
                full_text = "\n".join(param_lines) + "\n\n" + "\n".join(qe_cell_lines) + "\n\n" + "\n".join(atom_lines)
                QGuiApplication.clipboard().setText(full_text)
                QMessageBox.information(self, "コピー完了", f"スーパーセルの格子情報と原子座標 ({len(atoms)}個) をコピーしました。")
            except (ValueError, AttributeError, TypeError) as e:
                QMessageBox.critical(self, "コピーエラー", f"格子パラメータの取得または計算中にエラーが発生しました。\nスーパーセルが正しく生成されているか、元の格子パラメータが全て入力されているか確認してください。\n\n詳細: {e}"); return
        else:
            lines = [a.xyz_line() for a in atoms]
            QGuiApplication.clipboard().setText("\n".join(lines))
            QMessageBox.information(self, "コピー", f"Table {self.name} の内容 ({len(lines)}原子) をクリップボードにコピーしました。")
    def paste_clip(self):
        atoms = parse_xyz_block(QGuiApplication.clipboard().text())
        if atoms:
            atoms.sort(key=lambda a: a.z, reverse=True); self.load_atoms(atoms); self.main_window_ref.refresh()
        else: QMessageBox.warning(self, "貼り付けエラー", "クリップボードのデータ形式が不正か、データが含まれていません。")

HTML_TEMPLATE_FOR_APP="""<!DOCTYPE html><html><head><meta charset='utf-8'><script src='https://unpkg.com/ngl@2/dist/ngl.js'></script></head>
<body style='margin:0; overflow:hidden;'><div id='viewport' style='width:100vw; height:100vh;'></div>
<script>
var stage=new NGL.Stage('viewport',{backgroundColor:'white'});
function clearByName(name){stage.compList.slice().forEach(c=>{if(c.name===name)stage.removeComponent(c);});}
function loadPDB(name,str,colorVal,op,showLabel,showUnitcell){
 clearByName(name); var blob=new Blob([str],{type:'text/plain'});
 stage.loadFile(blob,{ext:'pdb',name:name}).then(function(c){
  var repParams={opacity:op,radiusScale:0.8}; if(colorVal&&colorVal.trim()!==""){repParams.color=colorVal;}
  c.addRepresentation('ball+stick',repParams);
  if(showLabel){c.addRepresentation('label',{labelType:'serial',color:'navy',scale:0.7,zOffset:1.5});}
  if(showUnitcell){c.addRepresentation('unitcell',{color:'grey'});}
  stage.autoView(200);
 });
}
function setBackground(color){stage.setParameters({backgroundColor:color});}
</script></body></html>"""

class Viewer(QWebEngineView):
    def __init__(self):
        super().__init__()
        self.temp_file = NamedTemporaryFile(delete=False, suffix='.html', mode='w', encoding='utf-8')
        self.temp_file.write(HTML_TEMPLATE_FOR_APP); self.temp_file.close()
        self.load(QUrl.fromLocalFile(self.temp_file.name)); self._ready = False
        self.loadFinished.connect(lambda _: setattr(self, '_ready', True))
    def update(self,name:str,pdb:str,color_val:Optional[str],op:float,show_label:bool,show_unitcell:bool):
        if self._ready:
            js_color_arg = f"'{color_val}'" if color_val else "null"
            safe_pdb = pdb.replace('\\', '\\\\').replace('`', '\\`')
            js_code = f"loadPDB('{name}', `{safe_pdb}`, {js_color_arg}, {op}, {str(show_label).lower()}, {str(show_unitcell).lower()});"
            self.page().runJavaScript(js_code)
    def remove(self,name:str):
        if self._ready: self.page().runJavaScript(f"clearByName('{name}');")
    def set_background(self, color: str):
        if self._ready: self.page().runJavaScript(f"setBackground('{color}');")
    def closeEvent(self, event):
        import os; os.unlink(self.temp_file.name); super().closeEvent(event)


class MainWindow(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle('Supercell Slab Structure Builder'); self.resize(1800, 900)
        self.supercell_na = 1; self.supercell_nb = 1
        self.supercell_va_orig: Optional[Tuple[float, ...]]=None; self.supercell_vb_orig: Optional[Tuple[float, ...]]=None; self.supercell_vc_orig: Optional[Tuple[float, ...]]=None
        self.tA = AtomTable('A', self); self.tB = AtomTable('B', self); self.tC = AtomTable('C (Supercell)', self); self.view = Viewer()
        main_layout = QHBoxLayout(self); left_panel = QVBoxLayout(); right_panel = QVBoxLayout()
        main_layout.addLayout(left_panel, 1); main_layout.addLayout(right_panel, 4)
        lattice_group = QGroupBox("Lattice Parameters"); lattice_layout = QGridLayout(); self.lattice_edits = {}
        params = ["a", "b", "c", "cosab", "cosac", "cosbc"]
        for i, param in enumerate(params):
            lattice_layout.addWidget(QLabel(f"{param}:"), i, 0); edit = QLineEdit(); self.lattice_edits[param] = edit; lattice_layout.addWidget(edit, i, 1)
        paste_params_btn = QPushButton("Paste Parameters"); paste_params_btn.clicked.connect(self.paste_lattice_params)
        paste_vectors_btn = QPushButton("Paste Vectors (CELL_PARAMETERS)"); paste_vectors_btn.clicked.connect(self.convert_vectors_to_params)
        lattice_layout.addWidget(paste_params_btn, len(params), 0, 1, 1); lattice_layout.addWidget(paste_vectors_btn, len(params), 1, 1, 1)
        lattice_group.setLayout(lattice_layout); left_panel.addWidget(lattice_group)
        supercell_group = QGroupBox("Supercell Builder"); supercell_layout = QGridLayout()
        supercell_layout.addWidget(QLabel("Sequence (e.g., AB,AA):"), 0, 0); self.sequence_edit = QLineEdit("AB,AA")
        supercell_layout.addWidget(self.sequence_edit, 0, 1); create_supercell_btn = QPushButton("Create Supercell")
        create_supercell_btn.clicked.connect(self.create_supercell); supercell_layout.addWidget(create_supercell_btn, 1, 0, 1, 2)
        supercell_group.setLayout(supercell_layout); left_panel.addWidget(supercell_group)
        display_group = QGroupBox("Display Options"); display_layout = QVBoxLayout()
        self.bg_combo = QComboBox(); self.bg_combo.addItems(["white", "black", "grey"]); self.bg_combo.currentTextChanged.connect(self.view.set_background)
        display_layout.addWidget(QLabel("Background Color:")); display_layout.addWidget(self.bg_combo)
        export_btn = QPushButton("3DビューをHTMLでエクスポート"); export_btn.clicked.connect(self.export_view_to_html)
        display_layout.addWidget(export_btn); display_group.setLayout(display_layout); left_panel.addWidget(display_group)
        self.controls = {}
        for name in ['A', 'B', 'C']:
            self.controls[name] = {'chk_visible':QCheckBox(f'表示 ({name})'),'chk_trans':QCheckBox(f'半透明 ({name})'),'chk_unitcell':QCheckBox(f'単位格子表示 ({name})')}
            self.controls[name]['chk_visible'].setChecked(True)
        for name in ['A', 'B', 'C']:
            for chk in self.controls[name].values():
                if isinstance(chk, QCheckBox): chk.stateChanged.connect(self.refresh)
        tables_layout = QHBoxLayout()
        tables_layout.addWidget(self._create_table_ui(self.tA, 'A', 'SkyBlue')); tables_layout.addWidget(self._create_table_ui(self.tB, 'B', 'Salmon'))
        tables_layout.addWidget(self._create_table_ui(self.tC, 'C', 'LightGreen')); right_panel.addLayout(tables_layout, 2)
        right_panel.addWidget(self.view, 3); left_panel.addStretch(1)

    def _create_table_ui(self, table_widget: AtomTable, name: str, default_color: str) -> QWidget:
        container=QWidget();layout=QVBoxLayout(container);layout.setContentsMargins(0,0,0,0)
        title=QLabel(f"<b>Table {table_widget.name}</b>");title.setAlignment(Qt.AlignCenter);layout.addWidget(title);layout.addWidget(table_widget)
        btn_layout_top=QHBoxLayout();paste_btn=QPushButton(f"Paste");copy_btn=QPushButton(f"Copy")
        paste_btn.clicked.connect(table_widget.paste_clip);copy_btn.clicked.connect(table_widget.copy_clip)
        btn_layout_top.addWidget(paste_btn);btn_layout_top.addWidget(copy_btn);layout.addLayout(btn_layout_top)
        if name in ['A','B']:
            fold_btn=QPushButton("単位格子内に折りたたむ");fold_btn.clicked.connect(lambda:self.fold_atoms_into_cell(table_widget));layout.addWidget(fold_btn)
        controls_frame=QFrame();controls_frame.setFrameShape(QFrame.StyledPanel);controls_layout=QVBoxLayout(controls_frame)
        controls_layout.addWidget(self.controls[name]['chk_visible']);controls_layout.addWidget(self.controls[name]['chk_trans'])
        controls_layout.addWidget(self.controls[name]['chk_unitcell']);controls_layout.addStretch();layout.addWidget(controls_frame)
        return container

    # ★★★★★ ここからが修正箇所 ★★★★★
    def paste_lattice_params(self, text=None):
        # text引数が文字列でない場合(ボタンクリックでFalseが渡された場合など)はクリップボードを使用する
        if not isinstance(text, str):
            text = QGuiApplication.clipboard().text()
        
        patterns={'a':r"a\s*=\s*(-?\d+\.?\d*e?[+-]?\d*)",'b':r"b\s*=\s*(-?\d+\.?\d*e?[+-]?\d*)",'c':r"c\s*=\s*(-?\d+\.?\d*e?[+-]?\d*)",
                  'cosab':r"cosab\s*=\s*(-?\d+\.?\d*e?[+-]?\d*)",'cosac':r"cosac\s*=\s*(-?\d+\.?\d*e?[+-]?\d*)",'cosbc':r"cosbc\s*=\s*(-?\d+\.?\d*e?[+-]?\d*)"}
        found_count=0
        for key,pattern in patterns.items():
            match=re.search(pattern,text,re.IGNORECASE)
            if match: self.lattice_edits[key].setText(match.group(1));found_count+=1
        if found_count>0:QMessageBox.information(self,"成功",f"{found_count}個の格子パラメータを貼り付けました。")
        else:QMessageBox.warning(self,"エラー","指定された書式の格子パラメータが見つかりませんでした。")
    # ★★★★★ ここまでが修正箇所 ★★★★★

    def convert_vectors_to_params(self):
        text=QGuiApplication.clipboard().text();lines=text.strip().splitlines();va,vb,vc=None,None,None
        try:
            header_index=-1
            for i,line in enumerate(lines):
                if"CELL_PARAMETERS"in line.upper():header_index=i;break
            if header_index==-1:raise ValueError("クリップボードに 'CELL_PARAMETERS' が見つかりません。")
            vector_lines=[ln for ln in lines[header_index+1:]if ln.strip()and not ln.strip().startswith('#')]
            if len(vector_lines)<3:raise ValueError("CELL_PARAMETERS の後に3行のベクトルデータが必要です。")
            va=tuple(map(float,vector_lines[0].split()));vb=tuple(map(float,vector_lines[1].split()));vc=tuple(map(float,vector_lines[2].split()))
            if any(len(v)!=3 for v in[va,vb,vc]):raise ValueError("各ベクトル行には3つの数値が必要です。")
        except(ValueError,IndexError)as e:QMessageBox.warning(self,"パースエラー",f"クリップボードの書式が不正です。\n{e}");return
        try:
            a=math.sqrt(sum(x*x for x in va));b=math.sqrt(sum(x*x for x in vb));c=math.sqrt(sum(x*x for x in vc))
            if a<1e-9 or b<1e-9 or c<1e-9:raise ValueError("格子ベクトルの長さがゼロです。")
            cosab=sum(x*y for x,y in zip(va,vb))/(a*b);cosac=sum(x*y for x,y in zip(va,vc))/(a*c);cosbc=sum(x*y for x,y in zip(vb,vc))/(b*c)
            self.lattice_edits['a'].setText(f"{a:.10f}");self.lattice_edits['b'].setText(f"{b:.10f}")
            self.lattice_edits['c'].setText(f"{c:.10f}");self.lattice_edits['cosab'].setText(f"{cosab:.10f}")
            self.lattice_edits['cosac'].setText(f"{cosac:.10f}");self.lattice_edits['cosbc'].setText(f"{cosbc:.10f}")
            QMessageBox.information(self,"変換成功","格子ベクトルからパラメータに変換し、値を入力しました。")
        except(ValueError,ZeroDivisionError)as e:QMessageBox.critical(self,"計算エラー",f"格子パラメータの計算中にエラーが発生しました。\n{e}");return
    def _calculate_lattice_vectors(self,params:Dict[str,float])->Optional[Tuple]:
        try:
            a,b,c=params['a'],params['b'],params['c'];cos_ab,cos_ac,cos_bc=params['cosab'],params['cosac'],params['cosbc']
            gamma=math.acos(cos_ab);beta=math.acos(cos_ac);alpha=math.acos(cos_bc);sin_gamma=math.sin(gamma)
            if abs(sin_gamma)<1e-9:raise ValueError("sin(gamma)がゼロに近すぎます。")
            va=(a,0.0,0.0);vb=(b*cos_ab,b*sin_gamma,0.0);cx=c*cos_ac;cy=c*(cos_bc-cos_ac*cos_ab)/sin_gamma
            vol_term=1-cos_bc**2-cos_ac**2-cos_ab**2+2*cos_bc*cos_ac*cos_ab
            if vol_term<-1e-9:raise ValueError("格子定数が物理的に不可能です。")
            cz=c*math.sqrt(max(0,vol_term))/sin_gamma;vc=(cx,cy,cz)
            return va,vb,vc
        except(ValueError,KeyError)as e:QMessageBox.critical(self,"格子ベクトル計算エラー",f"格子ベクトルの計算に失敗しました。\n{e}");return None
    def _invert_3x3_matrix(self,m:List[List[float]])->Optional[List[List[float]]]:
        det=(m[0][0]*(m[1][1]*m[2][2]-m[2][1]*m[1][2])-m[0][1]*(m[1][0]*m[2][2]-m[1][2]*m[2][0])+m[0][2]*(m[1][0]*m[2][1]-m[1][1]*m[2][0]))
        if abs(det)<1e-9:return None
        inv_det=1.0/det;inv_m=[[0,0,0]for _ in range(3)]
        inv_m[0][0]=(m[1][1]*m[2][2]-m[2][1]*m[1][2])*inv_det;inv_m[0][1]=(m[0][2]*m[2][1]-m[0][1]*m[2][2])*inv_det
        inv_m[0][2]=(m[0][1]*m[1][2]-m[0][2]*m[1][1])*inv_det;inv_m[1][0]=(m[1][2]*m[2][0]-m[1][0]*m[2][2])*inv_det
        inv_m[1][1]=(m[0][0]*m[2][2]-m[0][2]*m[2][0])*inv_det;inv_m[1][2]=(m[1][0]*m[0][2]-m[0][0]*m[1][2])*inv_det
        inv_m[2][0]=(m[1][0]*m[2][1]-m[2][0]*m[1][1])*inv_det;inv_m[2][1]=(m[2][0]*m[0][1]-m[0][0]*m[2][1])*inv_det
        inv_m[2][2]=(m[0][0]*m[1][1]-m[1][0]*m[0][1])*inv_det
        return inv_m
    def fold_atoms_into_cell(self,target_table:AtomTable):
        atoms=target_table.to_atoms()
        if not atoms:QMessageBox.information(self,"情報","テーブルに原子がありません。");return
        try:params={p:float(self.lattice_edits[p].text())for p in self.lattice_edits}
        except ValueError:QMessageBox.critical(self,"入力エラー","すべての格子パラメータが正しく入力されていません。");return
        vectors=self._calculate_lattice_vectors(params)
        if not vectors:return
        va,vb,vc=vectors
        m=[[va[i],vb[i],vc[i]]for i in range(3)];m_inv=self._invert_3x3_matrix(m)
        if not m_inv:QMessageBox.critical(self,"計算エラー","格子マトリクスが特異であり、逆行列を計算できません。");return
        folded_atoms=[]
        for atom in atoms:
            x,y,z=atom.x,atom.y,atom.z
            fx=m_inv[0][0]*x+m_inv[0][1]*y+m_inv[0][2]*z;fy=m_inv[1][0]*x+m_inv[1][1]*y+m_inv[1][2]*z;fz=m_inv[2][0]*x+m_inv[2][1]*y+m_inv[2][2]*z
            fx_new,fy_new,fz_new=fx%1.0,fy%1.0,fz%1.0
            x_new=m[0][0]*fx_new+m[0][1]*fy_new+m[0][2]*fz_new;y_new=m[1][0]*fx_new+m[1][1]*fy_new+m[1][2]*fz_new;z_new=m[2][0]*fx_new+m[2][1]*fy_new+m[2][2]*fz_new
            folded_atoms.append(Atom(atom.symbol,x_new,y_new,z_new,atom.f1,atom.f2,atom.f3))
        target_table.load_atoms(folded_atoms);self.refresh()
        QMessageBox.information(self,"成功",f"Table {target_table.name} の {len(folded_atoms)} 個の原子を単位格子内に折りたたみました。")
    def create_supercell(self):
        try:params={p:float(self.lattice_edits[p].text())for p in self.lattice_edits}
        except ValueError:QMessageBox.critical(self,"入力エラー","すべての格子パラメータが正しく入力されていません。");return
        vectors=self._calculate_lattice_vectors(params)
        if not vectors:return
        va,vb,vc=vectors;self.supercell_va_orig,self.supercell_vb_orig,self.supercell_vc_orig=va,vb,vc
        sequence_raw=self.sequence_edit.text().upper()
        if not sequence_raw.strip():QMessageBox.warning(self,"入力エラー","シーケンスが空です。");return
        b_axis_blocks=[block.strip()for block in sequence_raw.split(',')];self.supercell_nb=len(b_axis_blocks)if any(b_axis_blocks)else 1
        self.supercell_na=max([len(b)for b in b_axis_blocks]or[1]);atoms_A,atoms_B=self.tA.to_atoms(),self.tB.to_atoms()
        supercell_atoms:List[Atom]=[];b_axis_shift=(0.0,0.0,0.0)
        for block_str in b_axis_blocks:
            if not block_str:continue
            a_axis_shift=(0.0,0.0,0.0)
            for slab_type in block_str:
                source_atoms={'A':atoms_A,'B':atoms_B}.get(slab_type)
                if not source_atoms:QMessageBox.warning(self,"データエラー",f"シーケンス'{slab_type}'に対応するテーブルに原子データがありません。");return
                total_shift=tuple(x+y+z for x,y,z in zip(a_axis_shift,b_axis_shift,(0,0,0)))
                for atom in source_atoms:supercell_atoms.append(atom.get_shifted(total_shift))
                a_axis_shift=tuple(x+y for x,y in zip(a_axis_shift,va))
            b_axis_shift=tuple(x+y for x,y in zip(b_axis_shift,vb))
        if not supercell_atoms:QMessageBox.warning(self,"構築エラー","有効な原子が配置されませんでした。");return
        self.tC.load_atoms(supercell_atoms);QMessageBox.information(self,"成功",f"スーパーセルを構築しました ({len(supercell_atoms)} 原子)。");self.refresh()
    def refresh(self):
        tables={'A':self.tA,'B':self.tB,'C':self.tC};colors={'A':'skyblue','B':'salmon','C':'lightgreen'}
        for name,table in tables.items():
            controls=self.controls[name]
            if controls['chk_visible'].isChecked()and table.rowCount()>0:
                atoms=table.to_atoms();opacity=0.4 if controls['chk_trans'].isChecked()else 1.0
                show_unitcell=controls['chk_unitcell'].isChecked();lattice_params=None
                try:
                    if name in['A','B']:lattice_params={p:float(self.lattice_edits[p].text())for p in self.lattice_edits}
                    elif name=='C'and self.supercell_va_orig is not None:
                        a_o=float(self.lattice_edits['a'].text());b_o=float(self.lattice_edits['b'].text());c_o=float(self.lattice_edits['c'].text())
                        lattice_params={'a':a_o*self.supercell_na,'b':b_o*self.supercell_nb,'c':c_o,'cosab':float(self.lattice_edits['cosab'].text()),'cosac':float(self.lattice_edits['cosac'].text()),'cosbc':float(self.lattice_edits['cosbc'].text())}
                except(ValueError,KeyError):show_unitcell=False
                pdb_str=self._atoms_to_pdb(atoms,chain_id=name,lattice_params=lattice_params)
                self.view.update(name,pdb_str,colors[name],opacity,False,show_unitcell)
            else:self.view.remove(name)
    @staticmethod
    def _atoms_to_pdb(atoms:List[Atom],chain_id:str,lattice_params:Optional[Dict[str,float]]=None)->str:
        lines=[]
        if lattice_params:
            try:
                a,b,c=lattice_params['a'],lattice_params['b'],lattice_params['c']
                alpha=math.degrees(math.acos(lattice_params['cosbc']));beta=math.degrees(math.acos(lattice_params['cosac']));gamma=math.degrees(math.acos(lattice_params['cosab']))
                lines.append(f"CRYST1{a:9.3f}{b:9.3f}{c:9.3f}{alpha:7.2f}{beta:7.2f}{gamma:7.2f} P 1           1")
            except(KeyError,ValueError):pass
        lines.append(f"HEADER MODEL {chain_id} GENERATED BY SUPERCELL BUILDER")
        lines.extend(a.pdb_line(i+1,chain_id=chain_id)for i,a in enumerate(atoms));lines.append("END")
        return"\n".join(lines)
    def export_view_to_html(self):
        filePath,_=QFileDialog.getSaveFileName(self,"HTML形式でエクスポート","","HTML Files (*.html);;All Files (*)")
        if not filePath:return
        models_to_export=[]
        if self.controls['A']['chk_visible'].isChecked()and self.tA.rowCount()>0:models_to_export.append({'name':'A','table':self.tA,'color':'skyblue'})
        if self.controls['B']['chk_visible'].isChecked()and self.tB.rowCount()>0:models_to_export.append({'name':'B','table':self.tB,'color':'salmon'})
        if self.controls['C']['chk_visible'].isChecked()and self.tC.rowCount()>0:models_to_export.append({'name':'C','table':self.tC,'color':'lightgreen'})
        if not models_to_export:QMessageBox.information(self,"エクスポート","表示されている3Dモデルがありません。");return
        js_data_defs,js_promise_defs,js_promise_names=[],[],[]
        for model in models_to_export:
            name=model['name'];opacity=0.4 if self.controls[name]['chk_trans'].isChecked()else 1.0;show_uc=self.controls[name]['chk_unitcell'].isChecked();l_params=None
            try:
                if name in['A','B']:l_params={p:float(self.lattice_edits[p].text())for p in self.lattice_edits}
                elif name=='C'and self.supercell_va_orig is not None:
                    a_o,b_o,c_o=float(self.lattice_edits['a'].text()),float(self.lattice_edits['b'].text()),float(self.lattice_edits['c'].text())
                    l_params={'a':a_o*self.supercell_na,'b':b_o*self.supercell_nb,'c':c_o,'cosab':float(self.lattice_edits['cosab'].text()),'cosac':float(self.lattice_edits['cosac'].text()),'cosbc':float(self.lattice_edits['cosbc'].text())}
            except(ValueError,KeyError):show_uc=False
            pdb_str=self._atoms_to_pdb(model['table'].to_atoms(),chain_id=name,lattice_params=l_params).replace('`',r'\`')
            js_data_defs.append(f"const pdbData{name}=`{pdb_str}`;");js_data_defs.append(f"const blob{name}=new Blob([pdbData{name}],{{type:'text/plain'}});")
            js_promise_defs.append(f"""const promise{name}=stage.loadFile(blob{name},{{ext:'pdb',name:'Model{name}'}}).then(function(c){{
c.addRepresentation('ball+stick',{{color:'{model['color']}',opacity:{opacity},radiusScale:0.8}});
if({str(show_uc).lower()}){{c.addRepresentation('unitcell',{{color:'grey'}});}}}});""")
            js_promise_names.append(f"promise{name}")
        bg_color=self.bg_combo.currentText()
        script_content=f"""document.addEventListener('DOMContentLoaded',function(){{
var stage=new NGL.Stage('viewport',{{backgroundColor:'{bg_color}'}});
{' '.join(js_data_defs)} {' '.join(js_promise_defs)}
Promise.all([{','.join(js_promise_names)}]).then(function(){{stage.autoView();}});}});"""
        html_content=f"""<!DOCTYPE html><html lang="ja"><head><meta charset="utf-8"><title>3D Model Export</title>
<script src="https://unpkg.com/ngl@2/dist/ngl.js"></script></head><body style="margin:0; overflow:hidden;">
<div id="viewport" style="width:100vw; height:100vh;"></div><script>{script_content}</script></body></html>"""
        try:
            with open(filePath,'w',encoding='utf-8')as f:f.write(html_content)
            QMessageBox.information(self,"成功",f"HTMLファイルが正常に保存されました:\n{filePath}")
        except Exception as e:QMessageBox.critical(self,"保存エラー",f"ファイルの保存中にエラーが発生しました:\n{e}")
    def closeEvent(self,event):
        self.view.close();event.accept()

if __name__=='__main__':
    app=QApplication(sys.argv)
    win=MainWindow()
    win.show()
    sys.exit(app.exec_())