"""
    = vina.py =

    This plugin enables virtual screening with the AutoDock Vina software stack.
    It uses Meeko and Scrubber to prepare molecular ligands, and PLIP to
    analyze the results.
    
    It was tested on PyMOL 3.1 with Python 3.10. Currently supports only Linux
    and probably Mac.

    @author Pedro Sousa Lacerda
    @email pslacerda@gmail.com
"""


#
# SETUP CONDA PACKAGES
#
import subprocess
import shutil

try:
    import lxml, matplotlib, pandas, scipy, meeko, plip
except ImportError:
    subprocess.check_call(
        [
            "conda",
            "install",
            "-y",
            "-c",
            "conda-forge",
            "lxml",
            "matplotlib",
            "pandas",
            "openpyxl",
            "scipy",
            "meeko",
            "plip",
        ],
    )


#
# SETUP VINA
#

import platform
import subprocess
import sys
import os.path
from pymol import Qt

QStandardPaths = Qt.QtCore.QStandardPaths
data_dir = QStandardPaths.writableLocation(QStandardPaths.AppLocalDataLocation)
if not os.path.exists(data_dir):
    os.makedirs(data_dir)

system = platform.system().lower()
match system:
    case "windows":
        bin_fname = "vina_1.2.5_win.exe"
    case "linux":
        bin_fname = "vina_1.2.5_linux_x86_64"
    case "darwin":
        bin_fname = "vina_1.2.5_mac_x86_64"
vina_url = f"https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.5/{bin_fname}"
vina_bin = f"{data_dir}/vina"
if system == "windows":
    vina_bin += ".exe"

if not os.path.exists(vina_bin):
    import os
    import stat
    from urllib.request import urlretrieve
    print(f'Downloading "{vina_url}"')
    print(f'Installing Vina on "{vina_bin}"')
    urlretrieve(vina_url, vina_bin)
    os.chmod(vina_bin, stat.S_IEXEC)


#
# SETUP PIP
#

try:
    import scrubber
except ImportError:
    subprocess.check_call(
        [
            "python",
            "-m",
            "pip",
            "--disable-pip-version-check",
            "install",
            "https://github.com/forlilab/molscrub/archive/refs/tags/0.1.1.zip"
        ]
    )


#
# CODE STARTS HERE
#

import os
from os.path import (
    expanduser,
    dirname,
    splitext,
    basename,
)
from glob import glob
import itertools
from operator import itemgetter
import shutil
import shlex
import textwrap
import subprocess
import json
import atexit
import signal
from contextlib import contextmanager
import tempfile
from collections import Counter, OrderedDict
from pprint import pprint

import pymol
import pymol.gui
from pymol import cmd
from pymol.cgo import CYLINDER, SPHERE, COLOR
from pymol import Qt
import numpy as np
import pandas as pd
from lxml import etree
from matplotlib import pyplot as plt
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import euclidean

QWidget = Qt.QtWidgets.QWidget
QFileDialog = Qt.QtWidgets.QFileDialog
QFormLayout = Qt.QtWidgets.QFormLayout
QPushButton = Qt.QtWidgets.QPushButton
QSpinBox = Qt.QtWidgets.QSpinBox
QDoubleSpinBox = Qt.QtWidgets.QDoubleSpinBox
QDockWidget = Qt.QtWidgets.QDockWidget
QLineEdit = Qt.QtWidgets.QLineEdit
QCheckBox = Qt.QtWidgets.QCheckBox
QApplication = Qt.QtWidgets.QApplication
QMessageBox = Qt.QtWidgets.QMessageBox
QVBoxLayout = Qt.QtWidgets.QVBoxLayout
QTextEdit = Qt.QtWidgets.QTextEdit
QDialog = Qt.QtWidgets.QDialog
QDialogButtonBox = Qt.QtWidgets.QDialogButtonBox
QDesktopWidget = Qt.QtWidgets.QDesktopWidget
QProgressBar = Qt.QtWidgets.QProgressBar
QAction = Qt.QtWidgets.QAction
QComboBox = Qt.QtWidgets.QComboBox
QTabWidget = Qt.QtWidgets.QTabWidget
QTableWidget = Qt.QtWidgets.QTableWidget
QTableWidgetItem = Qt.QtWidgets.QTableWidgetItem
QHeaderView = Qt.QtWidgets.QHeaderView
QFrame = Qt.QtWidgets.QFrame
QDialogButtonBox = Qt.QtWidgets.QDialogButtonBox
QTreeView = Qt.QtWidgets.QTreeView

LeftDockWidgetArea = Qt.QtCore.Qt.LeftDockWidgetArea
QRegExp = Qt.QtCore.QRegExp
QtCore = Qt.QtCore
QThread = Qt.QtCore.QThread
pyqtSignal = Qt.QtCore.Signal
QStandardPaths = Qt.QtCore.QStandardPaths

QRegExpValidator = Qt.QtGui.QRegExpValidator
QPalette = Qt.QtGui.QPalette
QTextDocument = Qt.QtGui.QTextDocument
QIntValidator = Qt.QtGui.QIntValidator
QTextCursor = Qt.QtGui.QTextCursor
QIcon = Qt.QtGui.QIcon
QStandardItem = Qt.QtGui.QStandardItem
QStandardItemModel = Qt.QtGui.QStandardItemModel


###############################################
#                Utils                        #
###############################################

LIBRARIES_DIR = QStandardPaths.writableLocation(QStandardPaths.AppLocalDataLocation)
LIBRARIES_DIR += '/vina-pdbqt-libraries'
if not os.path.exists(LIBRARIES_DIR):
    os.makedirs(LIBRARIES_DIR)


TEMPDIR = tempfile.mkdtemp(prefix='pymol-vina-')

def clear_temp():
    shutil.rmtree(TEMPDIR)

atexit.register(clear_temp)

def run(command):
    ret = subprocess.run(
        shlex.split(command), stdout=subprocess.PIPE, stderr=subprocess.STDOUT
    )
    output = ret.stdout.decode()
    success = ret.returncode == 0
    return output, success


class BaseThread(QThread):

    numSteps = pyqtSignal(int)
    currentStep = pyqtSignal(int)

    logEvent = pyqtSignal(str)
    logCodeEvent = pyqtSignal(str)
    logRawEvent = pyqtSignal(str)

    done = pyqtSignal(bool)

    def __init__(self, *args, parent=None):
        super().__init__(parent)
        self.args = args


def display_box_sel(name, sel, margin):
    coords = cmd.get_coords(sel)
    max = np.max(coords, axis=0) + margin
    min = np.min(coords, axis=0) - margin
    display_box(name, max, min)


def display_box(name, max_coords, min_coords):
    #
    # From the original AutoDock plugin
    #

    box = [
        [max_coords[0], min_coords[0]],
        [max_coords[1], min_coords[1]],
        [max_coords[2], min_coords[2]],
    ]
    cylinder_size = 0.2
    color = [1.0, 1.0, 1.0]

    view = cmd.get_view()
    obj = []

    cmd.delete("_box")

    # box_color
    for i in range(2):
        for k in range(2):
            for j in range(2):
                if i != 1:
                    obj.append(CYLINDER)
                    obj.extend([box[0][i], box[1][j], box[2][k]])
                    obj.extend([box[0][i + 1], box[1][j], box[2][k]])
                    obj.append(cylinder_size)
                    obj.extend(color)
                    obj.extend(color)
                    obj.append(COLOR)
                    obj.extend(color)
                    obj.append(SPHERE)
                    obj.extend([box[0][i], box[1][j], box[2][k], cylinder_size])

                if j != 1:
                    obj.append(CYLINDER)
                    obj.extend([box[0][i], box[1][j], box[2][k]])
                    obj.extend([box[0][i], box[1][j + 1], box[2][k]])
                    obj.append(cylinder_size)
                    obj.extend(color)
                    obj.extend(color)
                    obj.append(COLOR)
                    obj.extend(color)
                    obj.append(SPHERE)
                    obj.extend([box[0][i], box[1][j + 1], box[2][k], cylinder_size])
                if k != 1:
                    obj.append(CYLINDER)
                    obj.extend([box[0][i], box[1][j], box[2][k]])
                    obj.extend([box[0][i], box[1][j], box[2][k + 1]])
                    obj.append(cylinder_size)
                    obj.extend(color)
                    obj.extend(color)
                    obj.append(COLOR)
                    obj.extend(color)
                    obj.append(SPHERE)
                    obj.extend([box[0][i], box[1][j], box[2][k + 1], cylinder_size])
    axes = [[2.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 2.0]]
    cmd.load_cgo(obj, name)
    cmd.set_view(view)



class OrderedCounter(Counter, OrderedDict):
    '''
    Counter that remembers the order elements are first encountered
    https://stackoverflow.com/a/23747652/199332
    '''
    def __repr__(self):
        return '%s(%r)' % (self.__class__.__name__, OrderedDict(self))

    def __reduce__(self):
        return self.__class__, (OrderedDict(self),)



###############################################
#          Load Result Pannel                 #
###############################################

class TreeItem(QStandardItem):
    def __init__(self, obj):
        super().__init__(str(obj))
        self.setFlags(self.flags() & ~QtCore.Qt.ItemIsEditable)


def parse_vina_pdbqt(filename):
    name = basename(filename)
    name = name.rsplit('.', maxsplit=2)[0]
        
    with open(filename) as pdbqt_file:
        for line in pdbqt_file:
            if line.startswith("MODEL"):
                _, mode_txt = line.split()
                mode = int(mode_txt)
            elif line.startswith("REMARK VINA RESULT:"):
                parts = line.split()
                affinity = float(parts[3])
                yield {
                    "name": name,
                    "filename": filename,
                    "affinity": affinity,
                    "mode": mode
                }

def load_plip_pose(project_dir, name, mode):
    cmd.delete("*")
    filename = f'{project_dir}/output/{name}.out.pdbqt'

    cmd.load(filename, 'lig', multiplex=True, zoom=False)
    cmd.set_name(f'lig_{mode.zfill(4)}', 'lig')
    cmd.delete('lig_*')
    cmd.alter('lig', 'chain="Z"')
    cmd.alter('lig', 'resn="LIG"')
    cmd.alter('lig', 'resi=1')

    filename = project_dir + '/target.pdbqt'
    cmd.load(filename, 'prot')

    pdb_fname = f"{TEMPDIR}/prot_lig.pdb"
    pse_fname = f'{TEMPDIR}/PROT_LIG_PROTEIN_LIG_Z_1.pse'

    cmd.save(pdb_fname, selection='*')
    command = f"conda run --cwd {TEMPDIR} plip -f {pdb_fname} -y --nohydro -o {TEMPDIR}"
    output, success = run(command)
    if not success:
        raise Exception(f"Failed to run PLIP: {command}")
    cmd.load(pse_fname)



def load_plip_full(project_dir, max_load, max_mode, tree_model):
    results = itertools.chain.from_iterable(
        map(parse_vina_pdbqt, glob(f"{project_dir}/output/*.out.pdbqt"))
    )
    results = list(sorted(results, key=lambda r: r['affinity']))
    cmd.set('pdb_conect_all', 'off')
    cmd.delete('prot')
    fname = f"{project_dir}/target.pdb"
    cmd.load(fname, 'prot')
    cmd.alter('prot', "type='ATOM'")
    fnames = []
    xml_l = []
    interactions = []
    interactions_type = [
        "hydrophobic_interaction",
        "hydrogen_bond",
        "water_bridge",
        "salt_bridge",
        "pi_stack",
        "pi_cation_interaction",
        "halogen_bond",
        "metal_complex"
    ]
    count = 0
    for idx, pose in enumerate(results):
        if pose['mode'] > max_mode:
            continue
        name = pose["name"]
        mode = str(pose["mode"])
        in_fname = project_dir + f'/output/{name}.out.pdbqt'
        out_fname = TEMPDIR + f'/{name}_m{mode}.out.pdb'
        
        cmd.delete('lig')
        cmd.load(in_fname, 'lig', multiplex=True, zoom=False)
        cmd.set_name(f'lig_{mode.zfill(4)}', 'lig')
        cmd.delete('lig_*')
        cmd.alter('lig', 'chain="Z"')
        cmd.alter('lig', 'resn="LIG"')
        cmd.alter('lig', 'resi=1')
        cmd.alter('lig', "type='HETATM'")
        cmd.save(out_fname, selection='*')
        command = f"conda run --cwd {TEMPDIR} plip -f '{out_fname}' -Oqsx --nohydro"
        proc = subprocess.run(shlex.split(command), stdout=subprocess.PIPE)
        assert proc.returncode == 0
        output = proc.stdout.decode().strip()
        plip = etree.fromstring(output)
        for inter_type in interactions_type:
            restype = plip.xpath(f"//{inter_type}/restype/text()")
            resnr = map(int, plip.xpath(f"//{inter_type}/resnr/text()"))
            reschain = plip.xpath(f"//{inter_type}/reschain/text()")
            for inter in zip(restype, resnr, reschain):
                interactions.append([f'{name}_m{mode}', inter_type, *inter])
        count += 1
        if count >= max_load:
            break
    
    interactions = sorted(interactions, key=lambda i: (i[4], i[3], i[1]))
    residues_l = ['%s%s%s' % (i[2], i[3], i[4]) for i in interactions]
    interactions_l = [i[1] for i in interactions]
    names_l = [i[0] for i in interactions]

    fig, axs = plt.subplots(len(interactions_type), layout="constrained", sharex=True)
    for ax, interaction_type in zip(axs, interactions_type):
        count = {}
        for res in residues_l:
            count[res] = 0
        for res, inter_type in zip(residues_l, interactions_l):
            if inter_type == interaction_type:
                count[res] += 1
        ax.bar(count.keys(), count.values())
        ax.set_title(interaction_type)
    plt.xticks(rotation=45)
    plt.show()

    fig, ax = plt.subplots(layout="constrained")
    df = pd.DataFrame({
        'name': names_l,
        'residue': residues_l
    })
    residues_l = list({r: None for r in residues_l})
    labels = []
    mols = []
    prev_name = None
    
    for cur_name, cur_residues in df.groupby('name'):
        if cur_name != prev_name:
            prev_name = cur_name
            counter = OrderedCounter(residues_l)
            for res in counter:
                counter[res] = 0
            for res in cur_residues['residue']:
                counter[res] += 1
            mols.append([counter[r] for r in residues_l])
            mol_item = TreeItem(cur_name)
            tree_model.appendRow(mol_item)
            for res, count in counter.items():
                if count > 0:
                    resn = res[:3]
                    chain = res[-1:]
                    resi = res[3:-1]
                    mol_item.appendRow([
                        TreeItem(chain),
                        TreeItem(resi),
                        TreeItem(resn),
                    ])
            labels.append(cur_name)
    
    D = []
    for idx1, mol1 in enumerate(mols):
        for idx2, mol2 in enumerate(mols):
            if idx1 >= idx2:
                continue
            d = euclidean(mol1, mol2) + 0.05
            D.append(d)
    Z = linkage(D)
    dendrogram(
        Z,
        labels=labels,
        orientation='right',
        color_threshold=0,
        distance_sort=True,
        ax=ax
    )
    ax.set_xlim(0)
    plt.show()


class ResultsWidget(QWidget):

    class ResultsTableWidget(QTableWidget):
        def __init__(self, project_dir):
            super().__init__()
            self.project_dir = project_dir
            self.props = ["Name", "Mode", "Affinity"]

            self.setSelectionBehavior(QTableWidget.SelectRows)
            self.setSelectionMode(QTableWidget.SingleSelection)
            self.setColumnCount(3)
            self.setHorizontalHeaderLabels(self.props)
            header = self.horizontalHeader()
            for idx in range(len(self.props)):
                header.setSectionResizeMode(
                    idx, QHeaderView.ResizeMode.ResizeToContents
                )

            @self.itemClicked.connect
            def itemClicked(item):
                name = self.item(item.row(), 0).text()
                mode = self.item(item.row(), 1).text()
                load_plip_pose(project_dir, name, mode)
        
    class SortableItem(QTableWidgetItem):
        def __init__(self, obj):
            super().__init__(str(obj))
            self.setFlags(self.flags() & ~QtCore.Qt.ItemIsEditable)

        def __lt__(self, other):
            try:
                return float(self.text()) < float(other.text())
            except ValueError:
                return self.text() < other.text()

    def __init__(self, is_intensive, project_dir, max_load, max_mode):
        super().__init__()
        self.is_intensive = is_intensive
        self.project_dir = project_dir
        self.max_load = max_load
        self.max_mode = max_mode

        layout = QVBoxLayout()
        self.setLayout(layout)

        tab = QTabWidget()
        layout.addWidget(tab)

        tab1_widget = QWidget()
        tab1_layout = QVBoxLayout(tab1_widget)
        tab1_widget.setLayout(tab1_layout)
        tab.addTab(tab1_widget, "Affinity list")

        self.table_widget = self.ResultsTableWidget(project_dir)
        tab1_layout.addWidget(self.table_widget)
        
        export_btn = QPushButton(QIcon("save"), "Export Table")
        export_btn.clicked.connect(self.export)
        tab1_layout.addWidget(export_btn)

        tab2_widget = QWidget()
        tab2_layout = QVBoxLayout(tab2_widget)
        tab2_widget.setLayout(tab2_layout)
        tab.addTab(tab2_widget, "Residue tree")
        
        self.tree_model = QStandardItemModel()
        self.tree_model.setHorizontalHeaderLabels(["Molecule/Chain", "Resi", "Resn"])
        self.tree_widget = QTreeView()
        self.tree_widget.setModel(self.tree_model)
        tab2_layout.addWidget(self.tree_widget) 
        
    def showEvent(self, event):
        self.refresh()
        super().showEvent(event)

    def refresh(self):
        self.table_widget.setSortingEnabled(False)

        # remove old rows
        while self.table_widget.rowCount() > 0:
            self.table_widget.removeRow(0)

        # append new rows
        project_dir = self.project_dir
        results = itertools.chain.from_iterable(
            map(parse_vina_pdbqt, glob(f"{self.project_dir}/output/*.out.pdbqt"))
        )
        results = sorted(results, key=itemgetter("affinity"))
        count = 0 
        for idx, pose in enumerate(results):
            if pose['mode'] <= self.max_mode:
                self.appendRow(pose)
                count += 1
            if count >= self.max_load:
                break
        self.table_widget.setSortingEnabled(True)

        if self.is_intensive:
            load_plip_full(self.project_dir, self.max_load, self.max_mode, self.tree_model)



    def appendRow(self, pose):
        self.table_widget.insertRow(self.table_widget.rowCount())
        line = self.table_widget.rowCount() - 1

        self.table_widget.setItem(line, 0, self.SortableItem(pose['name']))
        self.table_widget.setItem(line, 1, self.SortableItem(pose['mode']))
        self.table_widget.setItem(line, 2, self.SortableItem(pose['affinity']))
        
    def export(self):
        fileDialog = QFileDialog()
        fileDialog.setNameFilter("Excel file (*.xlsx)")
        fileDialog.setViewMode(QFileDialog.Detail)
        fileDialog.setAcceptMode(QFileDialog.AcceptSave)
        fileDialog.setDefaultSuffix(".xlsx")

        if fileDialog.exec_():
            filename = fileDialog.selectedFiles()[0]
            ext = os.path.splitext(filename)[1]
            with pd.ExcelWriter(filename) as xlsx_writer:
                row_count = self.table_widget.rowCount()
                col_count = self.table_widget.columnCount()
                data = []
                for row in range(row_count):
                    row_data = []
                    for col in range(col_count):
                        item = self.table_widget.item(row, col)
                        row_data.append(item.text() if item else '')
                    data.append(row_data)
                title = basename(self.project_dir)
                df = pd.DataFrame(data, columns=['Name', 'Mode', 'Affinity'])
                df.to_excel(xlsx_writer, sheet_name=title, index=False)
                  

def new_load_results_widget():
    dockWidget = QDockWidget()
    dockWidget.setWindowTitle("Analyze Vina")

    widget = QWidget()
    layout = QFormLayout(widget)
    widget.setLayout(layout)
    dockWidget.setWidget(widget)

    #
    # Max number of total loaded poses
    #
    max_load_spin = QSpinBox(widget)
    max_load_spin.setRange(1, 99999999)
    max_load_spin.setValue(15)
    max_load_spin.setGroupSeparatorShown(True)

    #
    # Only the best poses of each ligand
    #
    max_mode_spin = QSpinBox(widget)
    max_mode_spin.setRange(1, 20)
    max_mode_spin.setValue(9)
    max_mode_spin.setGroupSeparatorShown(True)

    #
    # Plot interaction histogram
    #
    intensive_check = QCheckBox()
    intensive_check.setChecked(False)

    #
    # Choose output folder
    #
    show_table_button = QPushButton("Load docking...", widget)

    @show_table_button.clicked.connect
    def load_results():
        nonlocal results_widget
        docking_file = str(
            QFileDialog.getOpenFileName(
                show_table_button,
                "Docking file",
                expanduser("~"),
                "Docking file (docking.json)",
            )[0]
        )
        if not docking_file:
            return
        
        project_dir = dirname(docking_file)

        if results_widget is not None:
            results_widget.setParent(None)
        del results_widget
        results_widget = ResultsWidget(
            intensive_check.isChecked(),
            project_dir,
            max_load_spin.value(),
            max_mode_spin.value(),
        )
        layout.setWidget(5, QFormLayout.SpanningRole, results_widget)

    #
    # Results Table
    #
    results_widget = None
    
    #
    # Setup form
    #
    layout.addRow("Max load:", max_load_spin)
    layout.addRow("Max mode:", max_mode_spin)
    layout.addRow("Intensive:", intensive_check)
    layout.setWidget(4, QFormLayout.SpanningRole, show_table_button)
    widget.setLayout(layout)

    return dockWidget


###############################################
#          Run Docking Pannel                 #
###############################################


class VinaThreadDialog(QDialog):
    def __init__(self, *vina_args, parent=None):
        super().__init__(parent)
        self.vina = VinaThread(*vina_args)
        self.vina.done.connect(self._done)

        # Setup window
        self.setModal(True)
        self.resize(QDesktopWidget().availableGeometry(self).size() * 0.7)
        self.setWindowFlags(self.windowFlags() | QtCore.Qt.CustomizeWindowHint)
        self.setWindowFlags(self.windowFlags() & ~QtCore.Qt.WindowCloseButtonHint)

        self.layout = QVBoxLayout(self)

        # Setup progress bar
        self.progress = QProgressBar()
        self.layout.addWidget(self.progress)
        self.progress.setValue(0)
        self.vina.numSteps.connect(self.progress.setMaximum)
        self.vina.currentStep.connect(self.progress.setValue)

        # Rich text output
        self.text = QTextEdit(self)
        self.layout.addWidget(self.text)
        self.text.setReadOnly(True)
        self.vina.logEvent.connect(self._appendHtml)
        self.vina.logCodeEvent.connect(self._appendCodeHtml)

        # Ok / Cancel buttons
        self.button_box = QDialogButtonBox(
            QDialogButtonBox.Ok | QDialogButtonBox.Abort, QtCore.Qt.Horizontal, self
        )
        self.layout.addWidget(self.button_box)
        self.button_box.accepted.connect(self._done)
        self.button_box.rejected.connect(self._abort)
        self.button_box.button(QDialogButtonBox.Ok).setDisabled(True)

        # Start docking
        self.vina.start()

    def _appendHtml(self, html):
        self.text.moveCursor(QTextCursor.End)
        self.text.insertHtml(self._prepareHtml(html))

    def _appendCodeHtml(self, html):
        self.text.moveCursor(QTextCursor.End)
        self.text.insertHtml("<pre>" + self._prepareHtml(html) + "</pre>")

    def _abort(self):
        self.vina.terminate()
        self.done(QDialog.Rejected)

    def _done(self, success):
        ok_button = self.button_box.button(QDialogButtonBox.Ok)
        abort_button = self.button_box.button(QDialogButtonBox.Abort)

        ok_button.setDisabled(False)
        abort_button.setDisabled(True)

        @self.button_box.accepted.connect
        def _done():
            if success:
                self.accept()
            else:
                self.reject()

    @staticmethod
    def _prepareHtml(html):
        return textwrap.dedent(html)


#
# Run docking software
#


class VinaThread(BaseThread):
    def run(self):
        (
            project_dir,
            ligands_file,
            target_sel,
            flex_sel,
            box_sel,
            box_margin,
            allow_errors,
            ph,
            exhaustiveness,
            num_modes,
            min_rmsd,
            energy_range,
            cpu,
            seed,
            save_library_check,
            library,
            scoring_func,
        ) = self.args

        #
        # Check previous output
        #
        if os.listdir(project_dir):
            self.logEvent.emit(f"""
                <br/>
                <font color="red">
                    <b>The docking folder is not empty: '{project_dir}'</b>
                </font>
            """)

        #
        # Prepare target
        #
        target_pdb = f"{project_dir}/target.pdb"
        target_basename = f"{project_dir}/target"
        cmd.save(target_pdb, target_sel)
        command = (
            f"conda run mk_prepare_receptor.py --read_pdb '{target_pdb}' -o '{target_basename}' -p"
        )
        if allow_errors:
            command = f"{command} -a"
        if flex_sel != "":
            flex_residues = set()
            for atom in cmd.get_model(flex_sel).atom:
                flex_residues.add(f"{atom.chain}:{atom.resi}")
            flex_residues = ",".join(flex_residues)
            command = f"{command} -f {flex_residues}"
        self.logEvent.emit(f"""
            <br/>
            <br/><b>Preparing target.</b>
            <br/><b>Command:</b> {command}
            <br/>
        """)
        output, success = run(command)
        self.logCodeEvent.emit(output)
        if not success:
            self.done.emit(False)
            return

        #
        # Create library
        #
        ligands_pdbqt = project_dir + "/ligands_pdbqt"

        if library:
            library_dir = LIBRARIES_DIR + '/' + library
            try:
                if os.path.exists(ligands_pdbqt):
                    shutil.rmtree(ligands_pdbqt)
            except OSError:
                os.unlink(ligands_pdbqt)
            os.symlink(library_dir, ligands_pdbqt, target_is_directory=True)
            self.logEvent.emit(f"""
                <br/>
                <br/><b>Using stored library:</b> {library_dir}
            """)
        elif ligands_file:
            if os.path.exists(ligands_pdbqt):
                try:
                    shutil.rmtree(ligands_pdbqt)
                except OSError:
                    os.unlink(ligands_pdbqt)
            #
            # Scrubbe isomers
            #
            ligands_sdf = project_dir + "/ligands.sdf"
            command = (
                f"conda run scrub.py -o '{ligands_sdf}' --ph {ph} --cpu {cpu} '{ligands_file}'"
            )
            self.logEvent.emit(
                f"""
                    <br/>
                    <br/><b>Scrubbing ligands.</b>
                    <br/><b>Command:</b> {command}
                    <br/>
                """
            )
            output, success = run(command)
            self.logCodeEvent.emit(output)
            if not success:
                self.done.emit(False)
                return

            #
            # Convert into PDBQT
            #
            command = (
                f"conda run mk_prepare_ligand.py -i '{ligands_sdf}' --multimol_outdir '{ligands_pdbqt}'"
            )
            self.logEvent.emit(f"""
                <br/>
                <br/><b>Converting ligands to PDBQT.</b>
                <br/><b>Command:</b> {command}
                <br/>
            """)
            output, success = run(command)
            self.logCodeEvent.emit(output)
            if not success:
                self.done.emit(False)
                return
            
            if save_library_check:
                library_dir = os.path.splitext(basename(ligands_file))[0]
                library_dir = LIBRARIES_DIR + '/' + library_dir
                self.logEvent.emit(f"""
                    <br/>
                    <br/><b>Storing library at:</b> {library_dir}
                """)
                try:
                    shutil.rmtree(library_dir)
                except:
                    pass
                shutil.copytree(ligands_pdbqt, library_dir)
        
        #
        # The number of dockings to do
        #
        count = len(glob(f"{ligands_pdbqt}/*.pdbqt"))
        n_ligands = count
        self.numSteps.emit(n_ligands)

        #
        # Compute box variables
        #
        box_coords = cmd.get_coords(box_sel)

        max = np.max(box_coords, axis=0)
        min = np.min(box_coords, axis=0)

        half_size = (max - min) / 2
        center = min + half_size

        size_x, size_y, size_z = (half_size + box_margin) * 2
        center_x, center_y, center_z = center

        size_x, size_y, size_z = (
            round(float(size_x), 2),
            round(float(size_y), 2),
            round(float(size_z), 2),
        )

        center_x, center_y, center_z = (
            round(float(center_x), 2),
            round(float(center_y), 2),
            round(float(center_z), 2),
        )
        
        #
        # Create Vina results directory
        #
        output_dir = f"{project_dir}/output"
        try:
            os.mkdir(output_dir)
        except FileExistsError:
            pass

        #
        # Project data
        #
        project_file = project_dir + "/docking.json"
        project_data = {
            "scoring_func": scoring_func,
            "size_x": size_x,
            "size_y": size_y,
            "size_z": size_z,
            "center_x": center_x,
            "center_y": center_y,
            "center_z": center_z,
        }

        #
        # Prompt for user confirmation
        #

        base_command = (
            f"{vina_bin}"
            f" --scoring {scoring_func}"
            f" --center_x {center_x}"
            f" --center_y {center_y}"
            f" --center_z {center_z}"
            f" --size_x {size_x}"
            f" --size_y {size_y}"
            f" --size_z {size_z}"
            f" --cpu {cpu}"
            f" --seed {seed}"
            f" --exhaustiveness {exhaustiveness}"
            f" --num_modes {num_modes}"
            f" --min_rmsd {min_rmsd}"
            f" --energy_range {energy_range}"
        )
        self.logEvent.emit(
            f"""
            <br/>
            <b>Vina base command:</b> {base_command}
        """
        )
        project_file = f"{project_dir}/docking.json"
        with open(project_file, "w") as docking_file:
            json.dump(project_data, docking_file, indent=4)

        fail_count = 0
        for idx, ligand_pdbqt in enumerate(glob(f"{ligands_pdbqt}/*.pdbqt")):
            name, _ = splitext(basename(ligand_pdbqt))
            output_pdbqt = f"{output_dir}/{name}.out.pdbqt"
            if os.path.exists(output_pdbqt):
                self.currentStep.emit(idx + 1)
                continue

            command = base_command + (
                f" --ligand '{ligand_pdbqt}'"
                f" --out '{output_pdbqt}'"
            )
            if os.path.exists(f"{project_dir}/flexible.pdbqt"):
                rigid_pdbqt = f"{project_dir}/rigid.pdbqt"
                flex_pdbqt = f"{project_dir}/flexible.pdbqt"
                command += f' --receptor "{rigid_pdbqt}"' f' --flex "{flex_pdbqt}"'
            else:
                target_pdbqt = f"{project_dir}/target.pdbqt"
                command += f' --receptor "{target_pdbqt}"'

            output, success = run(command)
            self.currentStep.emit(idx + 1)
            if not success:
                fail_count += 1
                if fail_count <= 10:
                    self.logEvent.emit(
                        f"""
                        <br/>
                        <font color="red">
                            <b>Vina command failed:</b> {command}
                            <br/>
                            <pre>{output}</pre>
                        </font>
                    """
                    )
                elif fail_count == 11:
                    self.logEvent.emit(
                        f"""
                        <br/>
                        <font color="red">
                            <b>Too many errors. Omitting output.</b>
                        </font>
                    """
                    )

        done_ligands = len(glob(f"{output_dir}/*.out.pdbqt"))

        self.logEvent.emit("<br/><h2>Summary</h2>")
        summary = f"""
            <br/><b>Total expected runs:</b> {n_ligands}
            <br/><b>Total failures:</b> {fail_count}
            <br/><b>Total found PDBQT files:</b> {done_ligands}
        """
        if done_ligands < n_ligands or fail_count > 0:
            self.logEvent.emit(f"<font color='red'>{summary}</font>")
        else:
            self.logEvent.emit(f"{summary}")

        self.done.emit(True)


def new_run_docking_widget():
    dockWidget = QDockWidget()
    dockWidget.setWindowTitle("Run Vina")

    widget = QWidget()

    layout = QFormLayout(widget)
    widget.setLayout(layout)
    dockWidget.setWidget(widget)

    #
    # Scoring function
    #
    scoring_func_combo = QComboBox(widget)
    scoring_func_combo.addItems(["vina", "vinardo"])

    #
    # Receptor selection
    #
    target_sel = QLineEdit("", widget)

    @target_sel.textEdited.connect
    def validate(text):
        validate_target_sel()

    def validate_target_sel():
        text = target_sel.text()
        palette = QApplication.palette(target_sel)
        palette.setColor(QPalette.Base, QtCore.Qt.white)
        valid = True
        try:
            if cmd.count_atoms(f"{text}") == 0:
                raise
        except:
            palette.setColor(QPalette.Base, QtCore.Qt.red)
            valid = False
        target_sel.setPalette(palette)
        return valid

    #
    # Flexible residues selection
    #
    flex_sel = QLineEdit("", widget)

    @flex_sel.textEdited.connect
    def validate(text):
        validate_flex_sel()

    def validate_flex_sel():
        text = flex_sel.text()
        palette = QApplication.palette(flex_sel)
        palette.setColor(QPalette.Base, QtCore.Qt.white)
        valid = True

        if text.strip() == "":
            palette.setColor(QPalette.Base, QtCore.Qt.white)
            return True
        try:
            if cmd.count_atoms(f"({text}) and ({target_sel.text()})") == 0:
                raise
            palette.setColor(QPalette.Base, QtCore.Qt.white)
            valid = True
        except:
            palette.setColor(QPalette.Base, QtCore.Qt.red)
            valid = False
        flex_sel.setPalette(palette)
        return valid

    #
    # Box selection
    #
    box_sel = QLineEdit("", widget)
    @box_sel.textEdited.connect
    def validate(text):
        validate_box_sel()
    def validate_box_sel():
        text = box_sel.text()
        palette = QApplication.palette(box_sel)
        palette.setColor(QPalette.Base, QtCore.Qt.white)
        try:
            if cmd.count_atoms(text) == 0:
                raise
        except:
            palette.setColor(QPalette.Base, QtCore.Qt.red)
            box_sel.setPalette(palette)
            cmd.delete("box")
            return False
        display_box_sel("box", text, box_margin_spin.value())
        box_sel.setPalette(palette)
        return True

    #
    # Miscellaneous options
    #
    box_margin_spin = QDoubleSpinBox(widget)
    box_margin_spin.setRange(0.0, 10.0)
    box_margin_spin.setValue(3.0)
    box_margin_spin.setSingleStep(0.1)
    box_margin_spin.setDecimals(1)
    @box_margin_spin.valueChanged.connect
    def display_box(margin):
        cmd.delete("box")
        display_box_sel("box", box_sel.text(), margin)

    allow_errors_check = QCheckBox(widget)
    allow_errors_check.setChecked(False)

    ph_spin = QDoubleSpinBox(widget)
    ph_spin.setRange(0.0, 14.0)
    ph_spin.setValue(7.0)
    ph_spin.setSingleStep(0.1)
    ph_spin.setDecimals(1)

    exhaustiveness_spin = QSpinBox(widget)
    exhaustiveness_spin.setRange(1, 50)
    exhaustiveness_spin.setValue(8)

    num_modes_spin = QSpinBox(widget)
    num_modes_spin.setRange(1, 20)
    num_modes_spin.setValue(9)

    min_rmsd_spin = QDoubleSpinBox(widget)
    min_rmsd_spin.setRange(0.0, 3.0)
    min_rmsd_spin.setValue(1.0)

    energy_range_spin = QDoubleSpinBox(widget)
    energy_range_spin.setRange(1.0, 10.0)
    energy_range_spin.setValue(3.0)

    cpu_count = QThread.idealThreadCount()
    cpu_spin = QSpinBox(widget)
    cpu_spin.setRange(1, cpu_count)
    cpu_spin.setValue(cpu_count)

    seed_spin = QSpinBox(widget)
    seed_spin.setRange(0, 10000)
    seed_spin.setValue(0)

    tab = QTabWidget()

    tab1_widget = QWidget()
    tab1_layout = QFormLayout(tab1_widget)
    tab1_widget.setLayout(tab1_layout)
    tab.addTab(tab1_widget, "New library")
    
    #
    # Choose ligand files and run docking
    #
    ligands_file = None
    ligands_button = QPushButton("Choose file...", widget)
    tab1_layout.addRow("Ligand file:", ligands_button)

    @ligands_button.clicked.connect
    def choose_ligands():
        nonlocal ligands_file
        ligands_file = str(
            QFileDialog.getOpenFileName(
                ligands_button, "Ligand files", expanduser("~"), "SMILES (*.smi);;SDF (*.sdf);;MOL (*.mol *.mol2)"
            )[0]
        )
        if not ligands_file:
            return
        ligands_button.setText(basename(ligands_file))
    #
    # Molecular library
    #
    save_library_check = QCheckBox()
    save_library_check.setChecked(False)
    tab1_layout.addRow("Store library:", save_library_check)


    tab2_widget = QWidget()
    tab2_layout = QFormLayout(tab2_widget)
    tab2_widget.setLayout(tab2_layout)
    tab.addTab(tab2_widget, "Stored library")

    tab_idx = 0
    library_combo = QComboBox()
    library_combo.addItems(os.listdir(LIBRARIES_DIR))
    tab2_layout.addRow("Library:", library_combo)
    @tab.currentChanged.connect
    def tab_changed(idx):
        nonlocal tab_idx
        tab_idx = idx
        if idx == 1:
            library_combo.clear()
            library_combo.addItems(os.listdir(LIBRARIES_DIR))

    #
    # Choose output folder
    #
    project_dir = None
    results_button = QPushButton("Choose folder...", widget)

    @results_button.clicked.connect
    def choose_project_dir():
        nonlocal project_dir
        project_dir = str(
            QFileDialog.getExistingDirectory(
                results_button,
                "Output folder",
                expanduser("~"),
                QFileDialog.ShowDirsOnly,
            )
        )
        if not project_dir:
            return
        results_button.setText(basename(project_dir))

    #
    # Run button
    #
    
    button = QPushButton("Run", widget)
    @button.clicked.connect
    def run():
        if not (validate_target_sel() & validate_flex_sel() & validate_box_sel()):
            return
        if not (project_dir):
            return
        nonlocal ligands_file
        library = library_combo.currentText().strip()
        if tab_idx == 0:
            if not ligands_file:
                return
            library = None
        elif tab_idx == 1:
            if not library:
                return
            ligands_file = None
        dialog = VinaThreadDialog(
            project_dir,
            ligands_file,
            target_sel.text(),
            flex_sel.text(),
            box_sel.text(),
            box_margin_spin.value(),
            allow_errors_check.isChecked(),
            ph_spin.value(),
            exhaustiveness_spin.value(),
            num_modes_spin.value(),
            min_rmsd_spin.value(),
            energy_range_spin.value(),
            cpu_spin.value(),
            seed_spin.value(),
            save_library_check.isChecked(),
            library,
            scoring_func_combo.currentText(),
        )
        dialog.exec_()

    horizontal_line = QFrame()
    horizontal_line.setFrameShape(QFrame.HLine)
    horizontal_line.setFrameShadow(QFrame.Sunken)

    #
    # setup layout
    #
    layout.addRow("Scoring function:", scoring_func_combo)
    layout.addRow("Target:", target_sel)
    layout.addRow("Flexible residues:", flex_sel)
    layout.addRow("Box:", box_sel)
    layout.addRow("Box margin:", box_margin_spin)
    layout.addRow("Allow Meeko errors:", allow_errors_check)
    layout.addRow("Ligand pH:", ph_spin)
    layout.addRow("Exhaustiveness:", exhaustiveness_spin)
    layout.addRow("Number of modes:", num_modes_spin)
    layout.addRow("Minimum RMSD:", min_rmsd_spin)
    layout.addRow("Energy range:", energy_range_spin)
    layout.addRow("Number of CPUs:", cpu_spin)
    layout.addRow("Seed number:", seed_spin)
    layout.setWidget(14, QFormLayout.SpanningRole, tab)
    layout.setWidget(15, QFormLayout.SpanningRole, horizontal_line)
    layout.addRow("Output folder:", results_button)
    layout.addWidget(button)
    widget.setLayout(layout)

    return dockWidget


def __init_plugin__(app=None):
    
    run_widget = new_run_docking_widget()
    load_widget = new_load_results_widget()

    run_widget.hide()
    load_widget.hide()

    window = pymol.gui.get_qtwindow()
    window.addDockWidget(LeftDockWidgetArea, run_widget)
    window.addDockWidget(LeftDockWidgetArea, load_widget)

    def show_run_widget():
        run_widget.show()
    
    def show_load_widget():
        load_widget.show()

    from pymol.plugins import addmenuitemqt
    addmenuitemqt("Vina (Run)", show_run_widget)
    addmenuitemqt("Vina (Analyze)", show_load_widget)


if __name__ in ["pymol", "pmg_tk.startup.XDrugPy"]:
    __init_plugin__()
