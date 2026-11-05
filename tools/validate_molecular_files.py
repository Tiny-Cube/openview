"""Acceptance checks against supplied molecules and real Gaussian text files."""
import argparse
from collections import Counter
import hashlib
import json
import multiprocessing
from pathlib import Path
import tempfile
import time
import numpy as np
from PySide6.QtCore import QTimer
from PySide6.QtTest import QTest
from PySide6.QtWidgets import QApplication,QMessageBox,QDoubleSpinBox
from openview.app import MainWindow,GaussianDialog,configure_graphics,STYLE
from openview.formats import load,gaussian_text
from openview.model import SYMBOLS

def file_hash(path):
    digest=hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):digest.update(block)
    return digest.hexdigest()

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--fchk',required=True);parser.add_argument('--chk',required=True);parser.add_argument('--log',required=True);parser.add_argument('--gjf',required=True);parser.add_argument('--out',default='analysis/molecular-acceptance.json');args=parser.parse_args()
    paths=[Path(args.fchk),Path(args.chk),Path(args.log),Path(args.gjf)];report={'files':[],'checks':[],'platform':__import__('sys').platform}
    hashes={p:file_hash(p) for p in paths}
    molecules=[]
    for path in paths:
        start=time.perf_counter();mol=load(path);molecules.append(mol)
        counts=Counter(SYMBOLS[int(z)] for z in mol.numbers)
        order=sorted(counts,key=lambda symbol:(0 if symbol=='C' else 1 if symbol=='H' else 2,symbol))
        formula=''.join(symbol+(str(counts[symbol]) if counts[symbol]>1 else '') for symbol in order)
        report['files'].append({'file':path.name,'sha256':hashes[path],'atoms':len(mol.numbers),'bonds':len(mol.bonds),'formula':formula,'charge':mol.charge,'multiplicity':mol.multiplicity,'energies':mol.energies,'frequencies':len(mol.frequencies),'frames':len(mol.frames),'bond_orders':np.unique(mol.bonds[:,2]).tolist(),'reader':mol.metadata.get('checkpoint_reader',path.suffix),'parse_seconds':time.perf_counter()-start})
    np.testing.assert_array_equal(molecules[0].numbers,molecules[1].numbers)
    np.testing.assert_allclose(molecules[0].frames,molecules[1].frames,atol=1e-6,rtol=0)
    np.testing.assert_array_equal(molecules[0].bonds,molecules[1].bonds)
    assert (molecules[0].charge,molecules[0].multiplicity)==(molecules[1].charge,molecules[1].multiplicity)
    report['checkpoint_coordinate_max_error_angstrom']=float(np.max(np.abs(molecules[0].frames-molecules[1].frames)))
    report['checks'].append('supplied CHK/FCHK atoms, positions, charge/spin and explicit bond orders agree')
    # Prove native CHK loading in isolation, with no companion FCHK and no formchk.
    with tempfile.TemporaryDirectory() as folder:
        isolated=Path(folder)/'standalone.chk';isolated.write_bytes(paths[1].read_bytes())
        direct=load(isolated);assert direct.metadata['checkpoint_reader']=='native-64'
        np.testing.assert_array_equal(direct.frames,molecules[1].frames)
    report['checks'].append('standalone native CHK, without companion FCHK')
    configure_graphics();app=QApplication([]);app.setStyle('Fusion');app.setStyleSheet(STYLE)
    window=MainWindow();window.show();QTest.qWait(200);errors=[]
    QMessageBox.warning=lambda parent,title,text:errors.append(title+': '+text)
    ticks=[];timer=QTimer();timer.setInterval(10);timer.timeout.connect(lambda:ticks.append(time.perf_counter()));timer.start()
    def wait(fn):
        end=time.perf_counter()+20
        while not fn() and time.perf_counter()<end:QTest.qWait(20)
        assert fn(),'Acceptance task timed out'
    output=Path(args.out);output.parent.mkdir(parents=True,exist_ok=True)
    for i,(path,reference) in enumerate(zip(paths,molecules)):
        start=time.perf_counter();before=len(ticks);window.open_path(path)
        wait(lambda:not window.loading and not window.tasks);QTest.qWait(80)
        assert not errors,errors
        mol=window.molecule;np.testing.assert_array_equal(mol.numbers,reference.numbers)
        assert window.view.ctx is not None and not window.view.error and window.view.draw_calls>=1
        assert window.table_model.rowCount()==len(mol.numbers)
        row=report['files'][i];row['gui_load_seconds']=time.perf_counter()-start;row['ui_heartbeats_during_load']=len(ticks)-before;row['gpu_draw_calls']=window.view.draw_calls
        if i==0:
            window.grab().save(str(output.with_suffix('.fchk.png')))
            window.tabs.setCurrentIndex(1);window.toggle_vibration();QTest.qWait(150)
            assert window.animation.isActive() and not np.array_equal(window.view.positions,mol.frames[-1])
            window.toggle_vibration();window.grab().save(str(output.with_suffix('.vibration.png')));window.tabs.setCurrentIndex(0)
        if i==1:window.grab().save(str(output.with_suffix('.chk.png')))
        window.pick(0);window.pick(1);assert 'Distance' in window.measurement.text()
        original_positions=mol.frames.copy();original_bonds=mol.bonds.copy();index=window.table_model.index(0,2)
        assert window.table_model.setData(index,str(float(mol.frames[window.frame,0,0])+.15))
        QTest.qWait(40);assert window.modified
        np.testing.assert_array_equal(mol.bonds,original_bonds)
        with tempfile.TemporaryDirectory() as folder:
            edited=Path(folder)/'edited.gjf';edited.write_text(gaussian_text(mol,window.frame,include_connectivity=True),encoding='utf-8');reloaded=load(edited)
            np.testing.assert_allclose(reloaded.frames[0],mol.frames[window.frame],atol=1e-6,rtol=0)
            # GJF groups connectivity by source atom; row order is not chemical identity.
            assert sorted(map(tuple,reloaded.bonds))==sorted(map(tuple,mol.bonds))
        window.undo();np.testing.assert_array_equal(mol.frames,original_positions);np.testing.assert_array_equal(mol.bonds,original_bonds)
        window.selection=[0];window.delete_atoms();QTest.qWait(20)
        assert len(mol.numbers)==len(reference.numbers)-1
        # Remaining explicit orders must survive index remapping after deletion.
        expected=original_bonds[(original_bonds[:,0]!=0)&(original_bonds[:,1]!=0)].copy();expected[:,:2]-=1
        np.testing.assert_array_equal(mol.bonds,expected)
        window.undo();np.testing.assert_array_equal(mol.bonds,original_bonds)
        window.modified=False
        row['edit_checks']=['table coordinate edit','original bond orders retained','GJF save/reopen','coordinate undo','atom delete/remap','topology undo']
    # Drive the actual Add atom dialog on the supplied molecule, then undo.
    window.accept_molecule(molecules[0]);original=window.molecule.frames.copy();bonds=window.molecule.bonds.copy()
    def fill_add_dialog():
        dialog=QApplication.activeModalWidget();assert dialog is not None
        fields=dialog.findChildren(QDoubleSpinBox)
        for field,value in zip(fields,[20.,20.,20.]):field.setValue(value)
        dialog.accept()
    QTimer.singleShot(40,fill_add_dialog);window.add_atom_coordinates();QTest.qWait(60)
    assert len(window.molecule.numbers)==len(original[0])+1
    np.testing.assert_array_equal(window.molecule.bonds,bonds);window.undo();np.testing.assert_array_equal(window.molecule.frames,original)
    report['checks'].extend(['real GPU display for FCHK/CHK/LOG/GJF','background parsing keeps Qt event loop active','selection and measurement','edit/export/reopen/undo for all four formats','add atom dialog and undo','FCHK vibration animation'])
    for path in paths:assert hashes[path]==file_hash(path)
    report['checks'].append('all original input files unchanged')
    report['gpu']=window.view.gpu_info;report['errors']=errors;timer.stop();window.modified=False;window.close()
    output.write_text(json.dumps(report,indent=2),encoding='utf-8');print(json.dumps(report,ensure_ascii=True))

if __name__=='__main__':
    multiprocessing.freeze_support();main()
