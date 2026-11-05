"""Additional native GUI checks for spectra, vibrations, surfaces and cancellation."""
import json
import multiprocessing
from pathlib import Path
import tempfile
import time
import numpy as np
from PySide6.QtWidgets import QApplication,QMessageBox
from PySide6.QtTest import QTest
from openview.app import MainWindow,GaussianDialog,configure_graphics,STYLE
from openview.model import Molecule
from openview.formats import gaussian_text

def main():
    configure_graphics();app=QApplication([]);app.setStyle('Fusion');app.setStyleSheet(STYLE)
    window=MainWindow();window.show();QTest.qWait(300);messages=[];checks=[]
    QMessageBox.warning=lambda parent,title,text:messages.append(text)
    def wait_for(fn,seconds=12):
        deadline=time.monotonic()+seconds
        while not fn() and time.monotonic()<deadline:QTest.qWait(20)
        assert fn(),'GUI task timeout'
    with tempfile.TemporaryDirectory() as folder:
        folder=Path(folder)
        coords=np.array([[0,0,0],[.9572,0,0],[-.24,.927,0]],dtype='f4')
        mol=Molecule([8,1,1],[coords,coords+.2],'Water · feature validation',frequencies=[-25.,1600.,3500.],intensities=[1.,30.,10.],modes=np.ones((3,3,3),dtype='f4')*.1,energies=[-76.,-76.1]).infer_bonds()
        window.accept_molecule(mol);window.tabs.setCurrentIndex(1);QTest.qWait(100)
        assert len(window.ir_plot.x)==2400 and len(window.energy_plot.y)==2
        window.toggle_vibration();QTest.qWait(160);assert window.mode_active and window.animation.isActive()
        assert not np.array_equal(window.view.positions,mol.frames[-1])
        window.toggle_vibration();assert not window.animation.isActive()
        window.toggle_animation();QTest.qWait(140);assert window.animation.isActive();window.toggle_animation()
        checks.extend(['IR broadening and plot','SCF energy plot','vibration displacement','trajectory playback'])
        dialog=GaussianDialog(mol,0,window);assert '%nprocshared=4' in dialog.preview.toPlainText();dialog.charge.setValue(1);assert not dialog.save.isEnabled();dialog.spin.setValue(2);assert dialog.save.isEnabled();dialog.close()
        checks.append('Gaussian input preview and electron parity')
        grid=np.indices((20,20,20)).astype('f4');values=(4.-np.sqrt(np.sum((grid-9.5)**2,axis=0))).astype('f4')
        mol.grid={'values':values,'axes':np.eye(3)*.2,'origin':np.array([-1.9]*3)}
        window.accept_molecule(mol);window.tabs.setCurrentIndex(2);window.iso.setValue(.3);window.build_surface()
        wait_for(lambda:not window.tasks)
        assert len(window.view.surfaces)==2;QTest.qWait(120);assert window.view.draw_calls==4
        out=Path('analysis');out.mkdir(exist_ok=True);window.grab().save(str(out/'gui-surfaces.png'))
        checks.extend(['background marching cubes','signed isosurface GPU rendering'])
        window.clear_surfaces();QTest.qWait(40);assert window.view.draw_calls==2
        # A superseded worker must not replace the most recently selected molecule.
        large=folder/'large.xyz';count=100000
        with large.open('w') as f:
            f.write(f'{count}\nCancel me\n')
            for i in range(count):f.write(f'C {i*2} 0 0\n')
        window.open_path(large);QTest.qWait(30);window.open_path(Path(__file__).resolve().parents[1]/'examples'/'water.xyz')
        wait_for(lambda:not window.loading and not window.tasks)
        assert len(window.molecule.numbers)==3
        checks.append('superseded parser cancellation')
        bad=folder/'bad.xyz';bad.write_text('2\nTruncated\nH 0 0 0\n');previous=window.molecule;window.open_path(bad)
        wait_for(lambda:not window.tasks);assert window.molecule is previous and messages
        checks.append('failed load retains current molecule')
        window.selection=[0,1];window.bond_order.setCurrentIndex(1);window.set_bond();QTest.qWait(50)
        assert window.view.bond_count==3
        window.bond_order.setCurrentIndex(2);window.set_bond();QTest.qWait(50);assert window.view.bond_count==4
        window.bond_order.setCurrentIndex(3);window.set_bond();QTest.qWait(50);assert window.view.bond_count==1
        window.undo();QTest.qWait(50);assert window.view.bond_count==4
        checks.append('double/triple/remove bond GPU instances and undo')
        window.modified=False;window.close()
        report={'checks':checks,'expected_invalid_file_message':messages,'unexpected_errors':[],'platform':__import__('sys').platform}
        (out/'gui-features.json').write_text(json.dumps(report,indent=2),encoding='utf-8');print(json.dumps(report))

if __name__=='__main__':
    multiprocessing.freeze_support();main()
