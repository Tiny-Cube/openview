import argparse
import copy
import multiprocessing
import os
import signal
from pathlib import Path
import sys
import time
import numpy as np
from PySide6.QtCore import Qt, QTimer, QProcess, QItemSelectionModel, QItemSelection
from PySide6.QtGui import QAction, QKeySequence, QSurfaceFormat, QFont
from PySide6.QtWidgets import (QApplication,QMainWindow,QWidget,QVBoxLayout,QHBoxLayout,
    QLabel,QPushButton,QFileDialog,QMessageBox,QComboBox,QSplitter,QTabWidget,QTableView,
    QSlider,QSpinBox,QDoubleSpinBox,QFormLayout,QDialog,QDialogButtonBox,QLineEdit,
    QPlainTextEdit,QGroupBox,QScrollArea,QHeaderView,QAbstractItemView,QCheckBox,QSizePolicy)
from .formats import save_xyz,gaussian_text
from .model import Molecule,SYMBOLS,element,measure,connected_fragment
from .renderer import MolecularView
from .table import AtomTable
from .plots import Plot,spectrum
from .work import Task
from .geometry_dialog import GeometryDialog
from .fragments import FRAGMENTS, build_element, element_options, fragment_name
from .fragment_palette import ElementPalette
from . import __version__

from .theme import STYLE,apply_theme


def button(text,callback,layout):
    b=QPushButton(text); b.clicked.connect(callback); layout.addWidget(b); return b

class GaussianDialog(QDialog):
    def __init__(self,mol,frame,parent):
        super().__init__(parent); self.setWindowTitle('Gaussian calculation input'); self.resize(720,640)
        self.mol=copy.copy(mol); self.frame=frame; self.saved_path=None
        layout=QVBoxLayout(self); form=QFormLayout(); layout.addLayout(form)
        self.charge=QSpinBox(); self.charge.setRange(-100,100); self.charge.setValue(mol.charge)
        self.spin=QSpinBox(); self.spin.setRange(1,99); self.spin.setValue(mol.multiplicity)
        self.cores=QSpinBox(); self.cores.setRange(1,1024); self.cores.setValue(4)
        self.memory=QLineEdit('2GB'); self.chk=QLineEdit('molecule.chk')
        self.route=QPlainTextEdit(mol.metadata.get('route','#p B3LYP/6-31G(d) Opt Freq')); self.route.setMaximumHeight(70)
        self.connectivity=QCheckBox('Include explicit bonds / bond orders');self.connectivity.setChecked(True)
        for name,field in [('Charge',self.charge),('Multiplicity',self.spin),('CPU cores',self.cores),('Memory',self.memory),('Checkpoint',self.chk),('Route',self.route)]:form.addRow(name,field)
        form.addRow(self.connectivity)
        self.preview=QPlainTextEdit(); self.preview.setReadOnly(True); self.preview.setFont(QFont('monospace',10)); layout.addWidget(self.preview)
        self.error=QLabel(); self.error.setWordWrap(True); layout.addWidget(self.error)
        row=QHBoxLayout(); layout.addLayout(row)
        self.copy_button=button('Copy input',lambda:QApplication.clipboard().setText(self.text()),row)
        self.save=button('Save .gjf',self.save_input,row); button('Close',self.reject,row)
        for field in (self.charge,self.spin,self.cores):field.valueChanged.connect(self.refresh)
        for field in (self.memory,self.chk):field.textChanged.connect(self.refresh)
        self.route.textChanged.connect(self.refresh); self.refresh()
        self.connectivity.toggled.connect(self.refresh)

    def text(self):
        self.mol.charge=self.charge.value(); self.mol.multiplicity=self.spin.value()
        return gaussian_text(self.mol,self.frame,self.route.toPlainText(),self.memory.text(),self.cores.value(),self.chk.text(),include_connectivity=self.connectivity.isChecked())
    def refresh(self):
        try:self.preview.setPlainText(self.text()); self.error.clear(); self.save.setEnabled(True);self.copy_button.setEnabled(True)
        except ValueError as exc:self.preview.clear(); self.error.setText(str(exc)); self.save.setEnabled(False);self.copy_button.setEnabled(False)
    def save_input(self):
        path,_=QFileDialog.getSaveFileName(self,'Save Gaussian input','molecule.gjf','Gaussian input (*.gjf *.com)')
        if path:
            try:
                Path(path).write_text(self.text(),encoding='utf-8')
                self.saved_path=Path(path); self.accept()
            except (OSError,ValueError) as exc:QMessageBox.warning(self,'Cannot save',str(exc))

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__(); apply_theme(QApplication.instance()); self.setWindowTitle('OpenView · Molecular Workbench'); self.resize(1380,860)
        self.molecule=None; self.frame=0; self.selection=[]; self.history=[]; self.redo_history=[]; self.modified=False
        self.tasks=[]; self.generation=0; self.edit_generation=0; self.surface_generation=0
        self.loading=False; self.job=None; self.job_offset=0; self.job_log=None
        self.animation=QTimer(self); self.animation.setTimerType(Qt.TimerType.PreciseTimer); self.animation.setInterval(33); self.animation.timeout.connect(self.animate)
        self.phase=0.; self.mode_active=False; self.frame_time=time.perf_counter(); self.file_path=None
        self.setAcceptDrops(True)
        toolbar=self.addToolBar('Molecule'); toolbar.setMovable(False)
        brand=QLabel('OpenView  '); brand.setObjectName('brand'); toolbar.addWidget(brand)
        self.action(toolbar,'Open',self.open_file,QKeySequence.StandardKey.Open)
        self.action(toolbar,'New molecule',self.new_molecule,QKeySequence.StandardKey.New)
        self.action(toolbar,'Example',self.example)
        self.action(toolbar,'Save GJF',self.gaussian,QKeySequence.StandardKey.Save)
        self.action(toolbar,'Export XYZ',self.export_xyz)
        self.action(toolbar,'Gaussian input',self.gaussian)
        self.action(toolbar,'Save PNG',self.screenshot)
        toolbar.addSeparator(); self.style=QComboBox(); self.style.addItems(['Ball & stick','Space filling','Sticks']); toolbar.addWidget(self.style)
        self.action(toolbar,'Fit [Home]',lambda:self.view.fit())
        self.action(toolbar,'Undo',self.undo,QKeySequence.StandardKey.Undo)
        self.action(toolbar,'Redo',self.redo,QKeySequence.StandardKey.Redo)
        root=QWidget(); self.setCentralWidget(root); outer=QVBoxLayout(root); outer.setContentsMargins(12,12,12,8)
        self.title=QLabel('Molecular workbench'); self.title.setStyleSheet('font-size:19px;font-weight:600;'); outer.addWidget(self.title)
        self.title.setSizePolicy(QSizePolicy.Policy.Ignored,QSizePolicy.Policy.Preferred)
        split=QSplitter(); outer.addWidget(split,1)
        center=QWidget(); cl=QVBoxLayout(center); cl.setContentsMargins(0,0,0,0)
        edit_tools=QHBoxLayout(); cl.addLayout(edit_tools)
        self.select_tool=button('Select / rotate',lambda:self.set_edit_mode('select'),edit_tools);self.select_tool.setCheckable(True);self.select_tool.setChecked(True)
        self.selection_tool=QComboBox();self.selection_tool.addItems(['Click','Rectangle','Brush'])
        self.selection_tool.setToolTip('Selection tool. Hold C/R/B for a temporary tool; F-click selects a connected fragment.')
        edit_tools.addWidget(self.selection_tool)
        self.add_tool=button('Build / replace',self.add_atom,edit_tools);self.add_tool.setCheckable(True)
        edit_tools.addWidget(QLabel('Element'))
        self.add_element=QComboBox();self.add_element.addItems(SYMBOLS[1:]);self.add_element.setCurrentText('C');self.add_element.setMinimumWidth(75)
        self.add_element.setMaxVisibleItems(12);self.add_element.activated.connect(lambda _:self.set_edit_mode('add'));edit_tools.addWidget(self.add_element)
        self.connect_added=QCheckBox('Bond new atom to selection');self.connect_added.setChecked(False);edit_tools.addWidget(self.connect_added);edit_tools.addStretch()
        self.carbon_fragment='tet'
        self.element_fragment='tet';self.element_fragments={'C':'tet'}
        self.fragment_palette=ElementPalette();self.carbon_palette=self.fragment_palette;cl.addWidget(self.fragment_palette)
        self.carbon_palette.set_fragment(self.carbon_fragment)
        self.fragment_palette.chosen.connect(self.choose_element_fragment)
        self.add_element.currentTextChanged.connect(self.update_fragment_palette)
        self.update_fragment_palette()
        self.edit_hint=QLabel('Click an atom to select. Choose an element and Add atom to place atoms.');self.edit_hint.setObjectName('muted');self.edit_hint.setWordWrap(True);cl.addWidget(self.edit_hint)
        self.edit_hint.setFixedHeight(self.edit_hint.fontMetrics().lineSpacing()*4+4)
        self.view=MolecularView(); cl.addWidget(self.view,1)
        controls=QHBoxLayout(); cl.addLayout(controls)
        self.play=button('▶ Play',self.toggle_animation,controls)
        self.timeline=QSlider(Qt.Orientation.Horizontal); self.timeline.setRange(0,0); self.timeline.valueChanged.connect(self.show_frame); controls.addWidget(self.timeline,1)
        self.frame_label=QLabel('Frame —'); controls.addWidget(self.frame_label)
        self.fps=QSpinBox(); self.fps.setRange(1,60); self.fps.setValue(30); self.fps.setSuffix(' fps'); controls.addWidget(self.fps)
        self.fps.valueChanged.connect(lambda v:self.animation.setInterval(round(1000/v)))
        hint=QLabel('Left: rotate · Middle/Shift: pan · Right/Ctrl: roll/zoom · Hold C/R/B: click/box/brush · F-click: fragment · Home: fit'); hint.setObjectName('muted'); hint.setWordWrap(True); cl.addWidget(hint)
        split.addWidget(center)
        self.tabs=QTabWidget(); self.tabs.setMinimumWidth(350); split.addWidget(self.tabs); split.setSizes([910,430])
        self.build_structure(); self.build_results(); self.build_surfaces(); self.build_jobs()
        self.view.atomPicked.connect(self.pick); self.view.graphicsError.connect(lambda s:self.statusBar().showMessage(s))
        self.view.selectionRequested.connect(self.apply_selection)
        self.selection_tool.currentIndexChanged.connect(self.change_selection_tool)
        self.view.placementRequested.connect(self.place_atom);self.view.exitPlacement.connect(lambda:self.set_edit_mode('select'))
        self.view.replacementRequested.connect(self.replace_atom)
        self.view.selectAllRequested.connect(self.select_all)
        self.view.clearSelectionRequested.connect(self.clear_selection)
        self.view.deleteRequested.connect(self.delete_atoms)
        self.style.currentTextChanged.connect(self.view.set_style)
        self.status=QLabel('Ready'); self.statusBar().addWidget(self.status,1)
        self.status.setSizePolicy(QSizePolicy.Policy.Ignored,QSizePolicy.Policy.Preferred)
        self.perf=QLabel(); self.statusBar().addPermanentWidget(self.perf)
        self.stats_timer=QTimer(self); self.stats_timer.setInterval(1000); self.stats_timer.timeout.connect(self.refresh_stats); self.stats_timer.start()
        self.build_menus(toolbar)
        self.example()

    def action(self,toolbar,text,fn,shortcut=None):
        action=QAction(text,self); action.triggered.connect(fn)
        if shortcut:action.setShortcut(shortcut)
        toolbar.addAction(action); return action

    def build_menus(self,toolbar):
        actions={a.text():a for a in toolbar.actions() if a.text()}
        file_menu=self.menuBar().addMenu('&File')
        for name in ('New molecule','Open','Example') :file_menu.addAction(actions[name])
        file_menu.addSeparator()
        for name in ('Save GJF','Export XYZ','Save PNG'):file_menu.addAction(actions[name])
        file_menu.addSeparator();self.action(file_menu,'Quit',self.close,QKeySequence.StandardKey.Quit)
        edit_menu=self.menuBar().addMenu('&Edit')
        edit_menu.addAction(actions['Undo']);edit_menu.addAction(actions['Redo'])
        redo_shortcuts=QKeySequence.keyBindings(QKeySequence.StandardKey.Redo)
        if QKeySequence('Ctrl+Shift+Z') not in redo_shortcuts:redo_shortcuts.append(QKeySequence('Ctrl+Shift+Z'))
        actions['Redo'].setShortcuts(redo_shortcuts)
        edit_menu.addSeparator();self.action(edit_menu,'Select all atoms',self.select_all)
        self.action(edit_menu,'Clear selection',self.clear_selection)
        self.action(edit_menu,'Delete selected atoms',self.delete_atoms)
        self.action(edit_menu,'Edit distance / angle / dihedral…',self.edit_geometry)
        view_menu=self.menuBar().addMenu('&View')
        view_menu.addAction(actions['Fit [Home]'])
        for index,name in enumerate(('Click selection','Rectangle selection','Brush selection')):
            self.action(view_menu,name,lambda checked=False,i=index:self.selection_tool.setCurrentIndex(i))
        help_menu=self.menuBar().addMenu('&Help')
        self.action(help_menu,'Mouse and keyboard controls',self.show_controls)
        self.action(help_menu,'About OpenView',lambda:QMessageBox.about(self,'About OpenView',
            f'OpenView {__version__}\nIndependent molecular workbench for Linux.\nTiny-Cube · 1620836684@qq.com\n\nGJF / LOG / FCHK; CHK supports identified layouts or an external formchk.\nGaussian is not included.'))

    def show_controls(self):
        QMessageBox.information(self,'OpenView controls',
            'View: left drag rotates; middle or Shift+left pans; right or Ctrl+left rolls and zooms. Wheel zooms; Home fits.\n\n'
            'Select: click atoms in measurement order. Hold C for clicks, R for a rectangle, B for a brush. Rectangle replaces the selection; Ctrl+rectangle adds. Brush adds atoms. Hold F and click for a connected fragment. A/N select all/none.\n\n'
            'Build: choose an element and Fragment, then click empty space to insert, H/open valence to attach, or a heavy atom to replace. Hold C to select without replacing. Esc returns to selection.\n\n'
            'Edit: select 2–4 atoms for distance/angle/dihedral editing. Set selected bond changes the bond order. Delete removes selected atoms; Ctrl+Z undoes, Ctrl+Shift+Z redoes.\n\n'
            'Save: Ctrl+S previews a GJF with explicit connections. XYZ does not store bond orders. Geometry edits clear stale calculation results.')

    def change_selection_tool(self,index):
        self.set_edit_mode('select');self.view.set_selection_mode(('click','rectangle','brush')[index]);self.view.setFocus()

    def build_structure(self):
        scroll=QScrollArea();scroll.setWidgetResizable(True)
        panel=QWidget(); layout=QVBoxLayout(panel);scroll.setWidget(panel); self.tabs.addTab(scroll,'Structure')
        self.summary=QLabel(); self.summary.setWordWrap(True); layout.addWidget(self.summary)
        self.measurement=QLabel('Select 2–4 atoms to measure'); self.measurement.setWordWrap(True); self.measurement.setStyleSheet('color:#096a63;padding:8px;background:#e1f3ef;border-radius:5px;'); layout.addWidget(self.measurement)
        self.geometry_button=button('Edit distance / angle / dihedral…',self.edit_geometry,layout)
        self.geometry_button.setToolTip('Select 2, 3 or 4 atoms in order, then adjust with live preview.')
        self.table_model=AtomTable(self); self.table=QTableView(); self.table.setModel(self.table_model)
        self.table.setMinimumHeight(140)
        self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setSelectionMode(QAbstractItemView.SelectionMode.ExtendedSelection)
        self.table.setAlternatingRowColors(True); self.table.verticalHeader().hide()
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        self.table.setColumnWidth(0,35); layout.addWidget(self.table,1)
        self.table_model.changed.connect(self.edit_atom)
        self.table.selectionModel().selectionChanged.connect(self.table_selection)
        row=QHBoxLayout(); layout.addLayout(row)
        button('Add by coordinates…',self.add_atom_coordinates,row); button('Delete selected',self.delete_atoms,row)
        row=QHBoxLayout(); layout.addLayout(row)
        button('Clear selection',self.clear_selection,row); button('Rebuild bonds',self.rebuild_bonds,row)
        row=QHBoxLayout(); layout.addLayout(row)
        self.bond_order=QComboBox();self.bond_order.addItems(['Single','Double','Triple','Remove']);row.addWidget(self.bond_order)
        button('Set selected bond',self.set_bond,row)
        note=QLabel('Double click a cell to edit. Bonds use covalent-distance inference; inferred bond orders are single.'); note.setObjectName('muted'); note.setWordWrap(True); layout.addWidget(note)

    def build_results(self):
        scroll=QScrollArea(); scroll.setWidgetResizable(True); panel=QWidget(); scroll.setWidget(panel); layout=QVBoxLayout(panel); self.tabs.addTab(scroll,'Results')
        self.result_info=QLabel(); self.result_info.setWordWrap(True); layout.addWidget(self.result_info)
        self.energy_plot=Plot('SCF energy'); layout.addWidget(self.energy_plot)
        self.ir_plot=Plot('IR spectrum · Lorentzian broadening'); layout.addWidget(self.ir_plot)
        form=QFormLayout(); layout.addLayout(form)
        self.mode=QComboBox(); self.mode.currentIndexChanged.connect(self.change_mode); form.addRow('Vibration',self.mode)
        self.amplitude=QDoubleSpinBox(); self.amplitude.setRange(0,3); self.amplitude.setSingleStep(.1); self.amplitude.setValue(.6); form.addRow('Displacement scale',self.amplitude)
        self.vibration_button=button('Animate vibration',self.toggle_vibration,layout)
        self.spectral_width=QDoubleSpinBox(); self.spectral_width.setRange(1,200); self.spectral_width.setValue(12); self.spectral_width.setSuffix(' cm⁻¹'); self.spectral_width.valueChanged.connect(self.update_spectrum); form.addRow('IR width',self.spectral_width)
        button('Export result CSV',self.export_results,layout); layout.addStretch()

    def build_surfaces(self):
        panel=QWidget(); layout=QVBoxLayout(panel); self.tabs.addTab(panel,'Surfaces')
        info=QLabel('Open a single-field Gaussian CUBE file. Build positive (green) and negative (violet) isosurfaces in a background process.'); info.setWordWrap(True); layout.addWidget(info)
        self.grid_info=QLabel('No scalar grid'); self.grid_info.setWordWrap(True); layout.addWidget(self.grid_info)
        form=QFormLayout(); layout.addLayout(form); self.iso=QDoubleSpinBox(); self.iso.setDecimals(6); self.iso.setRange(.000001,10000); self.iso.setValue(.03); self.iso.setSingleStep(.01); form.addRow('Absolute isovalue',self.iso)
        self.surface_button=button('Build ± isosurfaces',self.build_surface,layout)
        button('Clear surfaces',self.clear_surfaces,layout); layout.addStretch()

    def build_jobs(self):
        panel=QWidget(); layout=QVBoxLayout(panel); self.tabs.addTab(panel,'Jobs')
        info=QLabel('Run your separately installed Gaussian executable. Choose an executable and an input file; the job runs without blocking the interface.'); info.setWordWrap(True); layout.addWidget(info)
        form=QFormLayout(); layout.addLayout(form)
        self.executable=QLineEdit('g16'); form.addRow('Executable',self.executable)
        button('Choose executable',self.choose_executable,layout)
        self.job_input=QLineEdit(); form.addRow('Input .gjf / .com',self.job_input)
        button('Choose input',self.choose_job_input,layout)
        row=QHBoxLayout(); layout.addLayout(row); self.run_button=button('Run calculation',self.run_job,row); self.stop_button=button('Stop job',self.stop_job,row); self.stop_button.setEnabled(False)
        self.job_state=QLabel('No running job'); self.job_state.setWordWrap(True); layout.addWidget(self.job_state)
        self.job_output=QPlainTextEdit(); self.job_output.setReadOnly(True); self.job_output.setMaximumBlockCount(1200); layout.addWidget(self.job_output,1)
        self.job_timer=QTimer(self); self.job_timer.setInterval(500); self.job_timer.timeout.connect(self.tail_job)

    def confirm_discard(self):
        return not self.modified or QMessageBox.question(self,'Unsaved coordinates','Discard unsaved coordinate edits?',QMessageBox.StandardButton.Discard|QMessageBox.StandardButton.Cancel)==QMessageBox.StandardButton.Discard

    def open_file(self):
        path,_=QFileDialog.getOpenFileName(self,'Open molecule','','Molecular files (*.fchk *.fch *.chk *.log *.out *.gjf *.com *.xyz *.pdb *.ent *.mol *.sdf *.cube *.cub);;All files (*)')
        if path:self.open_path(path)

    def start_task(self,kind,args,on_success):
        # Superseded tasks are interrupted; retain each QThread until it finishes.
        for task in self.tasks:task.cancel()
        task=Task(kind,*args,parent=self); self.tasks.append(task)
        task.completed.connect(on_success); task.failed.connect(lambda message:self.worker_error(message) if not task.isInterruptionRequested() else None)
        task.finished.connect(lambda:self.task_finished(task)); task.start(); return task

    def task_finished(self,task):
        if task in self.tasks:self.tasks.remove(task)
        task.deleteLater()

    def worker_error(self,message):
        self.loading=False; self.status.setText('Operation failed; current molecule retained'); self.surface_button.setEnabled(self.molecule is not None and self.molecule.grid is not None)
        QMessageBox.warning(self,'Cannot complete operation',message)

    def open_path(self,path):
        if not self.confirm_discard():return
        self.animation.stop(); self.play.setText('▶ Play'); self.generation+=1; generation=self.generation
        self.loading=True; self.status.setText('Reading and building bonds in background…')
        def ready(mol):
            if generation!=self.generation:return
            self.file_path=Path(path); self.accept_molecule(mol); self.status.setText(f'Opened {Path(path).name}')
        self.start_task('load',(str(path),),ready)

    def example(self):
        if not self.confirm_discard():return
        self.generation+=1
        for task in self.tasks:task.cancel()
        # Original demonstration geometry; no proprietary fragments are reused.
        positions=[[1.396*np.cos(i*np.pi/3),1.396*np.sin(i*np.pi/3),0] for i in range(6)]
        positions.extend([[2.479*np.cos(i*np.pi/3),2.479*np.sin(i*np.pi/3),0] for i in range(6)])
        mol=Molecule([6]*6+[1]*6,positions,'Benzene · C₆H₆').infer_bonds(); self.file_path=None; self.accept_molecule(mol)

    def accept_molecule(self,mol):
        self.animation.stop(); self.loading=False; self.play.setText('▶ Play'); self.mode_active=False; self.phase=0
        self.set_edit_mode('select')
        self.frame=0; self.selection=[]; self.molecule=mol; self.history=[]; self.redo_history=[]; self.modified=False
        self.edit_generation+=1; self.surface_generation+=1; self.title.setText(mol.title)
        self.table_model.set_molecule(mol); self.view.set_molecule(mol); self.timeline.setRange(0,len(mol.frames)-1); self.timeline.setValue(0)
        self.mode.blockSignals(True); self.mode.clear()
        for i,f in enumerate(mol.frequencies): self.mode.addItem(f'{i+1} · {f:.2f} cm⁻¹'+(' · imaginary' if f<0 else ''))
        self.mode.blockSignals(False); self.vibration_button.setEnabled(mol.modes is not None and len(mol.modes)>0)
        self.energy_plot.set_data(np.arange(1,len(mol.energies)+1),mol.energies,'SCF evaluation','Hartree')
        self.update_spectrum(); self.update_summary(); self.show_frame(0)
        state=mol.metadata.get('normal_termination')
        self.result_info.setText(('Gaussian terminated normally' if state else 'Gaussian normal termination not found') if state is not None else 'Results appear when present in the opened file.')
        self.grid_info.setText(('Grid '+ ' × '.join(map(str,mol.grid['values'].shape))+f"\nRange {mol.grid['values'].min():.6g} to {mol.grid['values'].max():.6g}") if mol.grid else 'No scalar grid')
        self.surface_button.setEnabled(mol.grid is not None)
        if mol.metadata.get('source','').lower().endswith(('.log','.out')):
            self.timeline.setValue(len(mol.frames)-1)

    def update_summary(self):
        mol=self.molecule
        if mol is None:return
        self.summary.setText(f'{len(mol.numbers):,} atoms · {len(mol.bonds):,} bonds · {len(mol.frames):,} frames\nCharge {mol.charge:+d} · Multiplicity {mol.multiplicity}')
        self.measurement.setText(measure(mol.frames[self.frame],self.selection))
        self.geometry_button.setEnabled(not self.loading and 2<=len(self.selection)<=4)

    def edit_geometry(self):
        if self.loading or len(self.selection) not in (2,3,4):return
        self.animation.stop();self.mode_active=False;self.play.setText('▶ Play');self.show_frame(self.frame)
        for task in self.tasks:task.cancel()
        self.edit_generation+=1;self.surface_generation+=1
        mol=self.molecule;positions=mol.frames[self.frame]
        try:dialog=GeometryDialog(positions,mol.bonds,self.selection,self)
        except ValueError as exc:
            QMessageBox.warning(self,'Cannot edit geometry',str(exc));return
        surfaces=self.view.surfaces;self.view.surfaces=[];self.view.surface_dirty=True
        dialog.previewChanged.connect(self.view.set_positions)
        accepted=dialog.exec()==QDialog.DialogCode.Accepted
        self.view.surfaces=surfaces;self.view.surface_dirty=True;self.view.set_positions(positions)
        if accepted and not np.array_equal(dialog.candidate,positions):
            self.remember();mol.frames[self.frame]=dialog.candidate;self.invalidate_results()
            self.show_frame(self.frame);self.status.setText('Geometry updated · Ctrl+Z to undo')

    def show_frame(self,index):
        if self.molecule is None:return
        self.frame=index; self.table_model.refresh(index)
        self.view.set_positions(self.molecule.frames[index]); self.frame_label.setText(f'Frame {index+1} / {len(self.molecule.frames)}'); self.update_summary()

    def toggle_animation(self):
        if self.animation.isActive():self.animation.stop(); self.play.setText('▶ Play'); self.mode_active=False; self.show_frame(self.frame)
        elif self.molecule is not None and len(self.molecule.frames)>1:
            self.mode_active=False; self.animation.start(); self.play.setText('⏸ Pause')

    def animate(self):
        if self.mode_active:
            self.phase+=2*np.pi*self.animation.interval()/1000
            index=self.mode.currentIndex()
            if self.molecule.modes is None or index<0 or index>=len(self.molecule.modes):self.animation.stop();return
            self.view.set_positions(self.molecule.frames[-1]+self.amplitude.value()*np.sin(self.phase)*self.molecule.modes[index])
        else:self.timeline.setValue((self.frame+1)%len(self.molecule.frames))

    def toggle_vibration(self):
        if self.molecule is None or self.molecule.modes is None:return
        if self.mode_active and self.animation.isActive():self.animation.stop(); self.mode_active=False; self.play.setText('▶ Play'); self.show_frame(self.frame)
        else:self.timeline.setValue(len(self.molecule.frames)-1); self.mode_active=True; self.phase=0; self.animation.start(); self.play.setText('⏸ Pause')
    def change_mode(self):self.phase=0
    def update_spectrum(self):
        if self.molecule is None:return
        x,y=spectrum(self.molecule.frequencies,self.molecule.intensities,self.spectral_width.value()); self.ir_plot.set_data(x,y,'Wavenumber / cm⁻¹','IR intensity')

    def pick(self,index):
        if self.loading:return
        if index<0:self.selection=[]
        elif self.view.fragment_select_held:
            indices=connected_fragment(self.molecule.bonds,len(self.molecule.numbers),index)
            selected=set(self.selection);fragment=set(indices)
            self.selection=([i for i in self.selection if i not in fragment] if all(i in selected for i in indices)
                            else self.selection+[i for i in indices if i not in selected])
        elif index in self.selection:self.selection.remove(index)
        else:self.selection.append(index)
        self.sync_selection()

    def apply_selection(self,indices,additive=False):
        if self.loading or self.molecule is None:return
        selected=set(self.selection) if additive else set()
        result=list(self.selection) if additive else []
        for index in indices:
            index=int(index)
            if 0<=index<len(self.molecule.numbers) and index not in selected:
                result.append(index);selected.add(index)
        if result!=self.selection:self.selection=result;self.sync_selection()
    def sync_selection(self):
        sm=self.table.selectionModel(); sm.blockSignals(True); sm.clearSelection()
        # Batch contiguous ranges: selecting a large molecule must not issue N Qt updates.
        selected=sorted(set(self.selection)); ranges=QItemSelection()
        if selected:
            start=previous=selected[0]
            for index in selected[1:]+[None]:
                if index is not None and index==previous+1:previous=index;continue
                ranges.select(self.table_model.index(start,0),self.table_model.index(previous,4))
                start=previous=index
            sm.select(ranges,QItemSelectionModel.SelectionFlag.Select)
        sm.blockSignals(False); self.view.set_selection(self.selection); self.update_summary()
    def table_selection(self):
        self.selection=[i.row() for i in self.table.selectionModel().selectedRows()]; self.view.set_selection(self.selection); self.update_summary()
    def clear_selection(self):
        if self.loading:return
        self.selection=[];self.sync_selection()
    def select_all(self):
        if self.loading or self.molecule is None:return
        self.selection=list(range(len(self.molecule.numbers)));self.sync_selection()

    def replace_atom(self,index):
        if self.loading:return
        self.apply_element_fragment(target=index)

    def geometry_snapshot(self):
        mol=self.molecule
        size=mol.numbers.nbytes+mol.frames.nbytes+mol.bonds.nbytes
        return (mol.numbers.copy(),mol.frames.copy(),mol.bonds.copy(),size)

    def remember(self):
        mol=self.molecule
        # Geometry-only undo budget, at most 64 MiB and 20 edits.
        size=mol.numbers.nbytes+mol.frames.nbytes+mol.bonds.nbytes
        while self.history and (sum(x[3] for x in self.history)+size>64*1024**2 or len(self.history)>=20):self.history.pop(0)
        if size<=64*1024**2:self.history.append(self.geometry_snapshot())
        self.redo_history=[]
        self.animation.stop(); self.mode_active=False; self.play.setText('▶ Play'); self.modified=True

    def invalidate_results(self):
        # Geometry edits invalidate calculated fields; never display stale modes/results.
        self.molecule.modes=None; self.molecule.frequencies=[]; self.molecule.intensities=[]; self.molecule.energies=[]; self.molecule.grid=None
        self.vibration_button.setEnabled(False); self.mode.clear(); self.result_info.setText('Geometry edited. Recalculate to obtain matching results.'); self.energy_plot.set_data([],[]); self.ir_plot.set_data([],[])
        self.clear_surfaces(); self.surface_button.setEnabled(False); self.grid_info.setText('Grid invalidated by geometry edit')
        self.molecule.metadata={key:value for key,value in self.molecule.metadata.items() if key in ('source','route','input_coordinates')}

    def edit_atom(self,row,column,value):
        if self.loading:return
        for task in self.tasks:task.cancel()
        self.edit_generation+=1
        self.remember()
        if column==1:self.molecule.numbers[row]=value
        else:self.molecule.frames[self.frame,row,column-2]=value
        self.invalidate_results();self.view.set_positions(self.molecule.frames[self.frame]);self.update_summary();self.status.setText('Atom edited; existing bond orders preserved')

    def add_atom(self):
        self.set_edit_mode('add')

    def update_fragment_palette(self):
        symbol=self.add_element.currentText();number=element(symbol)
        self.element_fragment=self.element_fragments.get(symbol,element_options(number)['default'])
        self.fragment_palette.set_element(number,self.element_fragment)
        bare=self.element_fragment=='atom'
        self.connect_added.setEnabled(bare);self.connect_added.setVisible(bare)
        self.connect_added.setToolTip('Coordination fragments attach by clicking a hydrogen or an orange open valence.')
        if hasattr(self,'view') and self.view.interaction_mode=='add':self.set_edit_mode('add')

    def choose_carbon_fragment(self,key):
        self.add_element.setCurrentText('C');self.choose_element_fragment(key)

    def choose_element_fragment(self,key):
        symbol=self.add_element.currentText()
        self.element_fragments[symbol]=key;self.element_fragment=key
        if symbol=='C':self.carbon_fragment=key
        self.update_fragment_palette();self.set_edit_mode('add')

    def new_molecule(self):
        if not self.confirm_discard():return
        number=element(self.add_element.currentText())
        numbers,frames,bonds,_=build_element([],np.empty((1,0,3)),[],self.element_fragment,atomic_number=number,position=[0,0,0])
        mol=Molecule(numbers,frames,fragment_name(number,self.element_fragment),bonds=bonds)
        self.generation+=1
        for task in self.tasks:task.cancel()
        self.file_path=None;self.accept_molecule(mol);self.modified=True;self.set_edit_mode('add')

    def apply_carbon_fragment(self,*,position=None,target=None):
        self.apply_element_fragment(position=position,target=target)

    def apply_element_fragment(self,*,position=None,target=None):
        if self.loading or self.molecule is None:return
        mol=self.molecule
        try:
            numbers,frames,bonds,central=build_element(mol.numbers,mol.frames,mol.bonds,self.element_fragment,
                atomic_number=element(self.add_element.currentText()),position=position,target=target,orientation=self.view.rotation().T)
        except ValueError as exc:
            self.status.setText(str(exc));return
        if np.array_equal(numbers,mol.numbers) and np.array_equal(frames,mol.frames) and np.array_equal(bonds,mol.bonds):return
        self.remember();mol.numbers=numbers;mol.frames=frames;mol.bonds=bonds
        self.after_topology();self.selection=[central];self.sync_selection()
        remaining=int(np.count_nonzero(numbers==0))
        message=fragment_name(element(self.add_element.currentText()),self.element_fragment)+' · Ctrl+Z to undo'
        if remaining:message+=f' · {remaining} open valences: click an orange endpoint to extend'
        self.status.setText(message)

    def set_edit_mode(self,mode):
        if mode=='add' and (self.loading or self.molecule is None):return
        self.view.set_interaction_mode(mode)
        self.select_tool.setChecked(mode=='select');self.add_tool.setChecked(mode=='add')
        if mode=='add':
            self.animation.stop();self.mode_active=False;self.play.setText('▶ Play');self.show_frame(self.frame)
            self.edit_hint.setText(f'{self.add_element.currentText()}: click empty space to place a fragment; click H or an orange endpoint to extend; click a heavy atom to replace its coordination. Hold C: select · Esc: finish.')
        else:self.edit_hint.setText('Click atoms in measurement order; R-drag selects a box, B-drag brushes, F-click selects a connected fragment. Home fits the view. Choose an element to build; Ctrl+Z undoes.')

    def place_atom(self,position):
        if self.molecule is None or self.loading:return
        xyz=np.asarray(position,dtype='f4')
        if xyz.shape!=(3,) or not np.isfinite(xyz).all():return
        mol=self.molecule
        if np.any(np.linalg.norm(mol.frames[self.frame]-xyz,axis=1)<.2):
            self.status.setText('Too close to an existing atom; choose another position.');return
        if self.element_fragment!='atom':
            self.apply_element_fragment(position=xyz);return
        anchor=self.selection[0] if self.connect_added.isChecked() and len(self.selection)==1 else None
        self.remember();new_index=len(mol.numbers)
        mol.numbers=np.append(mol.numbers,element(self.add_element.currentText())).astype('i4')
        # Adding topology to a trajectory copies this position to all frames.
        mol.frames=np.concatenate((mol.frames,np.tile(xyz,(len(mol.frames),1,1))),axis=1)
        if anchor is not None:mol.bonds=np.vstack((mol.bonds,[anchor,new_index,1])).astype('f4')
        self.after_topology();self.selection=[new_index];self.sync_selection()
        text=f'Added {self.add_element.currentText()} at clicked position'
        if anchor is not None:text+=f' · single bond to atom {anchor+1}'
        self.status.setText(text+' · Esc to finish · Ctrl+Z to undo')

    def add_atom_coordinates(self):
        if self.molecule is None or self.loading:return
        dialog=QDialog(self); dialog.setWindowTitle('Add atom'); form=QFormLayout(dialog)
        symbol=QComboBox(); symbol.addItems(SYMBOLS[1:]); symbol.setCurrentText('C');form.addRow('Element',symbol)
        fields=[]
        for label in ('X / Å','Y / Å','Z / Å'):
            field=QDoubleSpinBox(); field.setRange(-100000,100000); field.setDecimals(5);form.addRow(label,field);fields.append(field)
        buttons=QDialogButtonBox(QDialogButtonBox.StandardButton.Ok|QDialogButtonBox.StandardButton.Cancel);form.addRow(buttons);buttons.accepted.connect(dialog.accept);buttons.rejected.connect(dialog.reject)
        if dialog.exec()!=QDialog.DialogCode.Accepted:return
        self.remember(); mol=self.molecule; mol.numbers=np.append(mol.numbers,element(symbol.currentText()))
        xyz=np.tile(np.array([f.value() for f in fields],dtype='f4'),(len(mol.frames),1,1))
        mol.frames=np.concatenate((mol.frames,xyz),axis=1); self.after_topology()

    def delete_atoms(self):
        if self.molecule is None or not self.selection or self.loading:return
        if len(self.selection)>=len(self.molecule.numbers):QMessageBox.information(self,'Keep one atom','A molecule must contain at least one atom.');return
        self.remember(); keep=np.ones(len(self.molecule.numbers),dtype=bool);keep[self.selection]=False
        bonds=self.molecule.bonds.copy();a=bonds[:,0].astype('i4');b=bonds[:,1].astype('i4')
        bonds=bonds[keep[a]&keep[b]];mapping=np.cumsum(keep)-1
        bonds[:,:2]=mapping[bonds[:,:2].astype('i4')];self.molecule.bonds=bonds
        self.molecule.numbers=self.molecule.numbers[keep];self.molecule.frames=self.molecule.frames[:,keep,:];self.after_topology()

    def after_topology(self):
        for task in self.tasks:task.cancel()
        self.edit_generation+=1;self.selection=[];self.invalidate_results()
        self.table_model.set_molecule(self.molecule);self.table_model.frame=self.frame;self.view.set_selection([])
        self.view.set_positions(self.molecule.frames[self.frame]);self.update_summary();self.status.setText('Topology edited; use Set selected bond to connect a new atom')

    def set_bond(self):
        if self.molecule is None or self.loading or len(self.selection)!=2:return
        for task in self.tasks:task.cancel()
        self.edit_generation+=1;self.remember();mol=self.molecule;a,b=sorted(self.selection)
        existing=mol.bonds;mask=~(((existing[:,0]==a)&(existing[:,1]==b))|((existing[:,0]==b)&(existing[:,1]==a)))
        mol.bonds=existing[mask]
        order=self.bond_order.currentIndex()+1
        if order<=3:mol.bonds=np.vstack((mol.bonds,[a,b,order])).astype('f4')
        self.invalidate_results();self.view.set_positions(mol.frames[self.frame]);self.update_summary();self.status.setText('Connectivity edited')

    def rebuild_bonds(self):
        if self.molecule is not None and np.any(self.molecule.numbers==0):
            self.status.setText('Finish or delete open valences before rebuilding bonds; their selected bond orders must be preserved.');return
        if self.molecule is None or self.loading:return
        self.animation.stop();self.mode_active=False;self.play.setText('▶ Play')
        self.edit_generation+=1; generation=self.edit_generation; mol=self.molecule
        self.view.set_positions(mol.frames[self.frame]); self.status.setText('Rebuilding bonds in background…')
        def ready(bonds):
            if mol is not self.molecule or generation!=self.edit_generation:return
            self.remember();mol.bonds=bonds;self.invalidate_results();self.view.set_positions(mol.frames[self.frame]);self.update_summary();self.status.setText('Bonds rebuilt · Ctrl+Z to restore previous bond orders')
        # Clear old connectivity immediately after topology edits.
        self.start_task('bonds',(mol.numbers,mol.frames[self.frame]),ready)

    def undo(self):
        if not self.history or self.loading:return
        self.redo_history.append(self.geometry_snapshot())
        self.restore_geometry(self.history.pop())

    def redo(self):
        if not self.redo_history or self.loading:return
        self.history.append(self.geometry_snapshot())
        self.restore_geometry(self.redo_history.pop())

    def restore_geometry(self,snapshot):
        for task in self.tasks:task.cancel()
        self.animation.stop();self.mode_active=False;self.play.setText('▶ Play')
        self.edit_generation+=1;z,frames,bonds,_=snapshot;self.molecule.numbers=z;self.molecule.frames=frames;self.molecule.bonds=bonds
        self.invalidate_results();self.selection=[];self.view.set_selection([]);self.table_model.set_molecule(self.molecule);self.show_frame(self.frame);self.modified=True

    def export_xyz(self):
        if self.molecule is None:return
        path,_=QFileDialog.getSaveFileName(self,'Export XYZ trajectory','molecule.xyz','XYZ (*.xyz)')
        if path:
            try:save_xyz(self.molecule,path);self.status.setText(f'Exported {path} · XYZ does not preserve bond orders; use Save GJF for connectivity')
            except (OSError,ValueError) as exc:QMessageBox.warning(self,'Cannot save',str(exc))
    def gaussian(self):
        if self.molecule is None or self.loading:return
        dialog=GaussianDialog(self.molecule,self.frame,self)
        if dialog.exec()==QDialog.DialogCode.Accepted and dialog.saved_path:
            same_state=(self.molecule.charge,self.molecule.multiplicity)==(dialog.mol.charge,dialog.mol.multiplicity)
            complete=len(self.molecule.frames)==1 and dialog.connectivity.isChecked() and same_state
            if complete:self.modified=False
            self.status.setText(f'Saved {dialog.saved_path}'+('' if complete else ' · calculation input exported; current molecule settings or trajectory are not fully saved'))
    def screenshot(self):
        path,_=QFileDialog.getSaveFileName(self,'Save viewport PNG','molecule.png','PNG (*.png)')
        if path:
            if not self.view.grabFramebuffer().save(path):QMessageBox.warning(self,'Cannot save','Failed to write image')
    def export_results(self):
        if self.molecule is None:return
        path,_=QFileDialog.getSaveFileName(self,'Export calculation results','results.csv','CSV (*.csv)')
        if not path:return
        import csv
        try:
            with open(path,'w',encoding='utf-8',newline='') as f:
                writer=csv.writer(f);writer.writerow(['SCF evaluation','Energy / Hartree'])
                writer.writerows(enumerate(self.molecule.energies,1));writer.writerow([]);writer.writerow(['Mode','Frequency / cm^-1','IR intensity'])
                for i,freq in enumerate(self.molecule.frequencies):writer.writerow([i+1,freq,self.molecule.intensities[i] if i<len(self.molecule.intensities) else ''])
        except OSError as exc:QMessageBox.warning(self,'Cannot save',str(exc))

    def build_surface(self):
        mol=self.molecule
        if mol is None or mol.grid is None or self.loading:return
        self.surface_generation+=1; generation=self.surface_generation;self.surface_button.setEnabled(False);self.status.setText('Building isosurfaces in background…')
        def ready(meshes):
            if mol is not self.molecule or generation!=self.surface_generation:return
            self.view.surfaces=meshes;self.view.surface_dirty=True;self.view.update();self.surface_button.setEnabled(True);self.status.setText(f'{sum(len(m[0])//3 for m in meshes):,} surface triangles')
        self.start_task('surface',(mol,self.iso.value()),ready)
    def clear_surfaces(self):
        self.surface_generation+=1;self.view.surfaces=[];self.view.surface_dirty=True;self.view.update()
        if self.molecule is not None:self.surface_button.setEnabled(self.molecule.grid is not None)

    def choose_executable(self):
        path,_=QFileDialog.getOpenFileName(self,'Choose Gaussian executable')
        if path:self.executable.setText(path)
    def choose_job_input(self):
        path,_=QFileDialog.getOpenFileName(self,'Choose Gaussian input','','Gaussian (*.gjf *.com)')
        if path:self.job_input.setText(path)
    def run_job(self):
        if self.job is not None:return
        path=Path(self.job_input.text()).expanduser()
        if not path.is_file():QMessageBox.warning(self,'Input required','Choose an existing Gaussian input file.');return
        output,_=QFileDialog.getSaveFileName(self,'Calculation output',str(path.with_suffix('.log')),'Gaussian log (*.log)')
        if not output:return
        if Path(output).resolve()==path.resolve():QMessageBox.warning(self,'Invalid output','Input and output must be different files.');return
        job=QProcess(self);job.setProgram(self.executable.text().strip());job.setWorkingDirectory(str(path.resolve().parent));job.setStandardInputFile(str(path.resolve()));job.setStandardOutputFile(output)
        if sys.platform.startswith('linux'):
            job.setUnixProcessParameters(QProcess.UnixProcessFlag.CreateNewSession)
        job.readyReadStandardError.connect(lambda:self.job_output.appendPlainText(bytes(job.readAllStandardError()).decode(errors='replace')))
        job.errorOccurred.connect(lambda error:self.job_error(job,error));job.finished.connect(lambda code,status:self.job_finished(job,code,status))
        self.job=job;self.job_log=Path(output);self.job_offset=0;self.job_output.clear();self.job_state.setText('Starting Gaussian…');self.run_button.setEnabled(False);self.stop_button.setEnabled(True)
        job.started.connect(lambda:self.job_state.setText(f'Running · PID {job.processId()}\n{output}'))
        job.start();self.job_timer.start()
    def tail_job(self):
        if self.job_log is None or not self.job_log.exists():return
        try:
            with self.job_log.open('rb') as f:
                f.seek(self.job_offset);data=f.read(65536);self.job_offset=f.tell()
            if data:self.job_output.appendPlainText(data.decode(errors='replace'))
        except OSError:pass
    def job_error(self,job,error):
        if self.job is job and error==QProcess.ProcessError.FailedToStart:self.finish_job_ui('Failed to start executable. Check its path and environment.')
    def job_finished(self,job,code,status):
        self.tail_job()
        if self.job is job:self.finish_job_ui(f'Process exited · code {code} · {status.name}\n{self.job_log}')
        job.deleteLater()
    def finish_job_ui(self,text):
        self.job_state.setText(text);self.job_timer.stop();self.run_button.setEnabled(True);self.stop_button.setEnabled(False)
        old=self.job;self.job=None
        if old is not None and old.state()==QProcess.ProcessState.NotRunning:old.deleteLater()
    def stop_job(self):
        if self.job is not None:
            job=self.job;pid=job.processId()
            if sys.platform.startswith('linux') and pid:
                try:os.killpg(pid,signal.SIGTERM)
                except ProcessLookupError:pass
            else:job.terminate()
            def force_stop():
                if self.job is job and job.state()!=QProcess.ProcessState.NotRunning:
                    if sys.platform.startswith('linux') and pid:
                        try:os.killpg(pid,signal.SIGKILL)
                        except ProcessLookupError:pass
                    else:job.kill()
            QTimer.singleShot(3000,force_stop)

    def refresh_stats(self):
        times=self.view.paint_times
        if times:self.perf.setText(f'GPU batches {self.view.draw_calls} · CPU paint {np.median(times[-60:]):.1f} ms')
    def dragEnterEvent(self,event):
        if event.mimeData().hasUrls() and event.mimeData().urls()[0].isLocalFile():event.acceptProposedAction()
    def dropEvent(self,event):self.open_path(event.mimeData().urls()[0].toLocalFile())
    def closeEvent(self,event):
        if self.job is not None:
            QMessageBox.information(self,'Calculation running','Stop the running calculation before closing.');event.ignore();return
        if not self.confirm_discard():event.ignore();return
        for task in self.tasks:task.cancel()
        for task in self.tasks:task.wait(3000)
        if any(task.isRunning() for task in self.tasks):event.ignore();return
        self.animation.stop();event.accept()

def configure_graphics():
    fmt=QSurfaceFormat();fmt.setVersion(3,3);fmt.setProfile(QSurfaceFormat.OpenGLContextProfile.CoreProfile);fmt.setDepthBufferSize(24);fmt.setSwapInterval(1);QSurfaceFormat.setDefaultFormat(fmt)

def main():
    multiprocessing.freeze_support()
    parser=argparse.ArgumentParser(description='OpenView molecular workbench');parser.add_argument('file',nargs='?');parser.add_argument('--version',action='version',version='OpenView '+__version__);args=parser.parse_args()
    configure_graphics();app=QApplication(sys.argv);app.setApplicationName('OpenView');app.setOrganizationName('OpenView');app.setStyle('Fusion');app.setStyleSheet(STYLE)
    window=MainWindow();window.show()
    if args.file:QTimer.singleShot(0,lambda:window.open_path(args.file))
    return app.exec()

if __name__=='__main__':sys.exit(main())
