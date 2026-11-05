"""Real mouse/pixel regression checks for light theme, selection and placement."""
import json
import multiprocessing
from pathlib import Path
import numpy as np
from PySide6.QtCore import Qt,QPoint,QPointF,QEvent,QTimer
from PySide6.QtGui import QImage,QPalette,QMouseEvent
from PySide6.QtTest import QTest
from PySide6.QtWidgets import QApplication,QMessageBox
from openview.app import MainWindow,configure_graphics
from openview.model import Molecule

def pixels(image):
    image=image.convertToFormat(QImage.Format.Format_RGBA8888)
    return np.frombuffer(image.constBits(),dtype=np.uint8).reshape(image.height(),image.bytesPerLine())[:,:image.width()*4].reshape(image.height(),image.width(),4).copy()

def main():
    configure_graphics();app=QApplication([]);window=MainWindow();window.show();QTest.qWait(150)
    errors=[];QMessageBox.warning=lambda p,t,s:errors.append(s)
    out=Path('analysis');out.mkdir(exist_ok=True);checks=[];view=window.view
    # A horizontal, separated C-O pair exposes the visible bond away from spheres.
    mol=Molecule([6,8],[[-1.5,0,0],[1.5,0,0]],'Selection regression',bonds=[[0,1,1]])
    window.accept_molecule(mol);view.yaw=0.;view.pitch=0.;view.update();QTest.qWait(100)
    assert view.ctx is not None and not view.error,view.error
    before=pixels(view.grabFramebuffer());xy,_=view.project()
    view.makeCurrent();bond_before=view.bond_resources[1].read();view.doneCurrent()
    QTest.mouseClick(view,Qt.MouseButton.LeftButton,pos=QPoint(*np.rint(xy[0]).astype(int)));QTest.qWait(80)
    assert window.selection==[0]
    after=pixels(view.grabFramebuffer());view.makeCurrent();bond_after=view.bond_resources[1].read();view.doneCurrent()
    assert bond_before==bond_after,'Selecting an atom changed bond GPU data'
    scale=view.devicePixelRatioF();cx,cy=np.rint(xy[0]*scale).astype(int)
    assert np.max(np.abs(after[cy,cx,:3].astype(int)-before[cy,cx,:3].astype(int)))>30,'Atom was not highlighted'
    # No changed framebuffer pixels beyond the selected atom footprint.
    yy,xx=np.indices(before.shape[:2]);radius=(.76*.48*view.height()/(2*view.scale))*scale+3
    outside=(xx-cx)**2+(yy-cy)**2>radius**2
    assert np.max(np.abs(after[outside].astype(int)-before[outside].astype(int)))<=1,'Selection leaked outside the atom'
    # Carbon half of the bond must have one stable color, without a gradient to oxygen.
    samples=[]
    for fraction in (.2,.3,.4):
        x,y=np.rint((xy[0]*(1-fraction)+xy[1]*fraction)*scale).astype(int);samples.append(after[y,x,:3].astype(int))
    assert np.max(np.ptp(np.asarray(samples),axis=0))<=2,'Bond color interpolates across its first half'
    view.grabFramebuffer().save(str(out/'interaction-selection.png'))
    checks.extend(['selection changes atom pixels only','bond GPU bytes unchanged by selection','crisp half-bond element colors'])
    # Popup palette and real popup capture, including highlighted text.
    window.add_element.showPopup();QTest.qWait(80);popup=window.add_element.view()
    palette=popup.palette();assert palette.color(QPalette.ColorRole.Base).lightness()>220
    assert palette.color(QPalette.ColorRole.Text).lightness()<100
    assert palette.color(QPalette.ColorRole.HighlightedText).lightness()<120
    popup.window().grab().save(str(out/'interaction-element-popup.png'));window.add_element.hidePopup()
    checks.append('light dropdown background and dark text, including selected rows')
    # Rotate/pan/zoom before placement, so screen coordinates cannot accidentally be treated as world XY.
    view.yaw=.67;view.pitch=-.41;view.pan[:]=[.28,-.33];view.scale=4.3;view.update();QTest.qWait(40)
    window.add_element.setCurrentText('N');window.choose_element_fragment('atom');QTest.mouseClick(window.add_tool,Qt.MouseButton.LeftButton)
    window.connect_added.setChecked(True)
    assert view.interaction_mode=='add'
    original=mol.frames.copy();old_bonds=mol.bonds.copy();original_origin=view.origin.copy();old_scale=view.scale
    point=QPoint(round(view.width()*.7),round(view.height()*.22))
    anchor_depth=view.project()[1][0]
    expected=view.screen_to_world(QPointF(point),float(anchor_depth))
    QTest.mouseClick(view,Qt.MouseButton.LeftButton,pos=point);QTest.qWait(60)
    assert len(mol.numbers)==3 and mol.numbers[-1]==7
    np.testing.assert_allclose(mol.frames[0,-1],expected,atol=1e-6)
    np.testing.assert_allclose(view.project()[0][-1],[point.x(),point.y()],atol=.01)
    np.testing.assert_array_equal(mol.bonds,[[0,1,1],[0,2,1]])
    np.testing.assert_array_equal(view.origin,original_origin);assert view.scale==old_scale
    assert window.selection==[2] and window.table_model.rowCount()==3
    window.undo();QTest.qWait(30);np.testing.assert_array_equal(mol.frames,original);np.testing.assert_array_equal(mol.bonds,old_bonds)
    checks.extend(['clicked placement follows rotated/panned/zoomed view and anchor depth','optional single bond to selected anchor','new atom selected for chaining','camera preserved on add','undo restores geometry and bonds'])
    # Hold C to select without replacing in build mode, as in the original GV.
    xy,_=view.project();QTest.keyPress(view,Qt.Key.Key_C)
    QTest.mouseClick(view,Qt.MouseButton.LeftButton,pos=QPoint(*np.rint(xy[1]).astype(int)))
    QTest.keyRelease(view,Qt.Key.Key_C)
    assert len(mol.numbers)==2 and window.selection==[1]
    assert mol.numbers[1]==8
    QTest.mouseClick(view,Qt.MouseButton.LeftButton,pos=QPoint(*np.rint(xy[1]).astype(int)))
    assert len(mol.numbers)==2 and mol.numbers[1]==7
    window.undo();assert mol.numbers[1]==8
    # Dragging rotates but cannot place an atom. Explicit move event includes held button.
    start=QPoint(40,40);end=QPoint(120,100);yaw=view.yaw
    QTest.mousePress(view,Qt.MouseButton.LeftButton,pos=start)
    event=QMouseEvent(QEvent.Type.MouseMove,QPointF(end),QPointF(view.mapToGlobal(end)),Qt.MouseButton.NoButton,Qt.MouseButton.LeftButton,Qt.KeyboardModifier.NoModifier)
    app.sendEvent(view,event);QTest.mouseRelease(view,Qt.MouseButton.LeftButton,pos=end);QTest.qWait(30)
    assert len(mol.numbers)==2 and view.yaw!=yaw
    QTest.keyClick(view,Qt.Key.Key_Escape);assert view.interaction_mode=='select' and window.select_tool.isChecked()
    QTest.mouseClick(view,Qt.MouseButton.LeftButton,pos=QPoint(30,30));assert len(mol.numbers)==2 and not window.selection
    checks.extend(['C-click selects anchor without replacing','build click replaces element and undo restores it','drag rotates without adding','Esc returns to selection; empty click clears'])
    # Original GV camera mappings: right/Ctrl drag roll and zoom, middle/Shift pan.
    def drag(button,modifiers):
        QTest.mousePress(view,button,modifiers,pos=start)
        event=QMouseEvent(QEvent.Type.MouseMove,QPointF(end),QPointF(view.mapToGlobal(end)),Qt.MouseButton.NoButton,button,modifiers)
        app.sendEvent(view,event);QTest.mouseRelease(view,button,modifiers,pos=end)
    pan=view.pan.copy();roll=view.roll;scale=view.scale
    drag(Qt.MouseButton.RightButton,Qt.KeyboardModifier.NoModifier)
    assert view.roll!=roll and view.scale!=scale
    np.testing.assert_array_equal(view.pan,pan)
    roll=view.roll;scale=view.scale
    drag(Qt.MouseButton.LeftButton,Qt.KeyboardModifier.ControlModifier)
    assert view.roll!=roll and view.scale!=scale
    scale=view.scale;roll=view.roll
    drag(Qt.MouseButton.MiddleButton,Qt.KeyboardModifier.NoModifier)
    assert not np.array_equal(view.pan,pan) and view.scale==scale and view.roll==roll
    QTest.keyClick(view,Qt.Key.Key_A);assert window.selection==[0,1]
    QTest.keyClick(view,Qt.Key.Key_N);assert not window.selection
    checks.extend(['right and Ctrl drag roll/zoom; middle drag pans','A/N select all/none'])
    # With bonding off, placing an atom keeps old connections even when an anchor is selected.
    window.selection=[0];window.sync_selection();window.connect_added.setChecked(False);window.add_element.setCurrentText('O');window.choose_element_fragment('atom');window.add_atom()
    QTest.mouseClick(view,Qt.MouseButton.LeftButton,pos=QPoint(80,80));QTest.qWait(40)
    assert len(mol.numbers)==3;np.testing.assert_array_equal(mol.bonds,old_bonds);window.undo()
    checks.append('bonding-off placement preserves existing topology')
    # Exercise the real geometry dialog, including preview and one-step undo.
    window.accept_molecule(Molecule([6,8],[[-1.,0,0],[1.,0,0]],bonds=[[0,1,2]]))
    window.selection=[0,1];window.sync_selection();original=window.molecule.frames.copy()
    observed=[]
    def change_geometry():
        dialog=QApplication.activeModalWidget()
        dialog.value.setFocus();dialog.value.selectAll();QTest.keyClicks(dialog.value,'2.5');QTest.keyClick(dialog.value,Qt.Key.Key_Tab)
        QTest.qWait(40)
        observed.append(float(np.linalg.norm(view.positions[1]-view.positions[0])))
        dialog.grab().save(str(out/'geometry-dialog.png'))
        dialog.accept()
    QTimer.singleShot(50,change_geometry);QTest.mouseClick(window.geometry_button,Qt.MouseButton.LeftButton)
    assert abs(observed[0]-2.5)<1e-5 and len(window.history)==1
    assert abs(float(np.linalg.norm(window.molecule.frames[0,1]-window.molecule.frames[0,0]))-2.5)<1e-5
    window.undo();np.testing.assert_array_equal(window.molecule.frames,original)
    checks.append('geometry dialog keyboard edit, preview, commit and undo')
    assert not errors,errors
    window.modified=False;window.example();view.yaw=.3;view.pitch=-.25;view.update();window.connect_added.setChecked(True);window.add_element.setCurrentText('C');window.set_edit_mode('add');window.status.setText('Ready · click to place carbon');QTest.qWait(50)
    window.grab().save(str(out/'interaction-light-ui.png'))
    report={'checks':checks,'errors':errors,'platform':__import__('sys').platform,'gpu':view.gpu_info}
    (out/'interaction-check.json').write_text(json.dumps(report,indent=2),encoding='utf-8');print(json.dumps(report));window.modified=False;window.close()

if __name__=='__main__':multiprocessing.freeze_support();main()
