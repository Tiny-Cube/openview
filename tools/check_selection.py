"""Real Qt selection acceptance for disconnected structures and build mode."""
import argparse
import json
import multiprocessing
from pathlib import Path
import numpy as np
from PySide6.QtCore import Qt, QPoint, QPointF, QEvent
from PySide6.QtGui import QMouseEvent, QFocusEvent
from PySide6.QtTest import QTest
from PySide6.QtWidgets import QApplication, QMessageBox
from openview.app import MainWindow, configure_graphics
from openview.model import Molecule


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--out',default='validation/0.4.0/selection.json');args=parser.parse_args()
    target=Path(args.out);target.parent.mkdir(parents=True,exist_ok=True)
    configure_graphics();app=QApplication([]);window=MainWindow();window.show();QTest.qWait(150)
    errors=[];QMessageBox.warning=lambda p,t,s:errors.append(s)
    mol=Molecule([6,8,7,1,16,26],[[-3,0,0],[-1.5,0,0],[0,0,0],[2,0,0],[3,0,0],[0,3,0]],
                 'Selection acceptance',bonds=[[0,1,2],[1,2,1],[3,4,1]])
    window.accept_molecule(mol);view=window.view;view.yaw=view.pitch=view.roll=0.;view.origin[:]=0;view.scale=5.;view.pan[:]=0;view.update();QTest.qWait(80)
    assert view.ctx is not None and not view.error,view.error
    coordinates=mol.frames.copy();bonds=mol.bonds.copy();checks=[]
    def point(index):return QPoint(*np.rint(view.project()[0][index]).astype(int))
    def move(position,modifiers=Qt.KeyboardModifier.NoModifier):
        e=QMouseEvent(QEvent.Type.MouseMove,QPointF(position),QPointF(view.mapToGlobal(position)),Qt.MouseButton.NoButton,Qt.MouseButton.LeftButton,modifiers)
        app.sendEvent(view,e)
    def box(indices,add=False,release_key_first=False,temporary=True):
        xy=view.project()[0][indices]
        start=QPoint(*np.floor(xy.min(axis=0)-12).astype(int));end=QPoint(*np.ceil(xy.max(axis=0)+12).astype(int))
        modifiers=Qt.KeyboardModifier.ControlModifier if add else Qt.KeyboardModifier.NoModifier
        if temporary:QTest.keyPress(view,Qt.Key.Key_R,modifiers)
        QTest.mousePress(view,Qt.MouseButton.LeftButton,modifiers,pos=start);move(end,modifiers)
        assert view.rubber_band.isVisible()
        if temporary and not add:
            QTest.qWait(30);window.grab().save(str(target.with_name('selection-rectangle-held.png')))
        if temporary and release_key_first:QTest.keyRelease(view,Qt.Key.Key_R)
        QTest.mouseRelease(view,Qt.MouseButton.LeftButton,modifiers,pos=end)
        if temporary and not release_key_first:QTest.keyRelease(view,Qt.Key.Key_R)
        QTest.qWait(30);assert not view.rubber_band.isVisible() and view.selection_drag is None
    yaw=view.yaw;box([0,1,2],release_key_first=True)
    assert window.selection==[0,1,2] and view.yaw==yaw
    assert len(window.table.selectionModel().selectedRows())==3
    checks.append('R rectangle commits on mouse release even when R is released first;camera and table selection agree')
    scale=view.scale;box([3,4],add=True)
    assert window.selection==[0,1,2,3,4] and view.scale==scale
    checks.append('Ctrl+rectangle adds another disconnected group without zooming')
    window.clear_selection();QTest.keyPress(view,Qt.Key.Key_B)
    QTest.mousePress(view,Qt.MouseButton.LeftButton,pos=point(3));move(point(4));move(point(3))
    QTest.mouseRelease(view,Qt.MouseButton.LeftButton,pos=point(3));QTest.keyRelease(view,Qt.Key.Key_B)
    assert window.selection==[3,4] and view.yaw==yaw
    checks.append('B brush adds each encountered atom once without toggling or rotating')
    window.clear_selection();window.selection_tool.setCurrentIndex(2)
    QTest.keyPress(view,Qt.Key.Key_F);QTest.mouseClick(view,Qt.MouseButton.LeftButton,pos=point(0))
    assert window.selection==[0,1,2] and view.scale==scale
    QTest.mouseClick(view,Qt.MouseButton.LeftButton,pos=point(1));QTest.keyRelease(view,Qt.Key.Key_F)
    assert not window.selection
    window.selection_tool.setCurrentIndex(0)
    checks.append('F-click toggles the connected component and leaves disconnected atoms and camera unchanged')
    window.add_element.setCurrentText('N');window.choose_element_fragment('atom')
    QTest.keyPress(view,Qt.Key.Key_F);QTest.mouseClick(view,Qt.MouseButton.LeftButton,pos=point(0));QTest.keyRelease(view,Qt.Key.Key_F)
    assert window.selection==[0,1,2] and mol.numbers[0]==6
    QTest.keyClick(view,Qt.Key.Key_Escape)
    checks.append('F-click in build mode selects a fragment without replacing the clicked atom')
    window.clear_selection();QTest.keyPress(view,Qt.Key.Key_R)
    QTest.mousePress(view,Qt.MouseButton.LeftButton,pos=QPoint(10,10));move(QPoint(view.width()-10,view.height()-10))
    app.sendEvent(view,QFocusEvent(QEvent.Type.FocusOut))
    QTest.mouseRelease(view,Qt.MouseButton.LeftButton,pos=QPoint(view.width()-10,view.height()-10))
    assert not window.selection and not view.rubber_band.isVisible() and not view.held_selection_modes
    checks.append('focus loss cancels an uncommitted rectangle and clears temporary modes')
    window.selection_tool.setCurrentIndex(1);box([0,1,2],temporary=False)
    assert window.selection==[0,1,2]
    window.selection_tool.setCurrentIndex(0);view.pan[:]=[1.,1.];QTest.keyClick(view,Qt.Key.Key_Home)
    np.testing.assert_array_equal(view.pan,[0.,0.])
    assert len(window.menuBar().actions())==4
    checks.append('persistent rectangle tool,Home fit and File/Edit/View/Help menus available')
    np.testing.assert_array_equal(mol.frames,coordinates);np.testing.assert_array_equal(mol.bonds,bonds)
    assert not window.modified and not errors,errors
    checks.append('selection never changes geometry,connectivity or the dirty flag')
    # QAction shortcuts require the owning window to be active, unlike direct widget key handlers.
    window.activateWindow();window.raise_();QApplication.setActiveWindow(window);view.setFocus();QTest.qWait(50)
    assert window.isActiveWindow() and view.hasFocus()
    window.edit_atom(0,2,float(coordinates[0,0,0])+.2)
    QTest.keyClick(view,Qt.Key.Key_Z,Qt.KeyboardModifier.ControlModifier);QTest.qWait(30)
    np.testing.assert_array_equal(mol.frames,coordinates)
    QTest.keyClick(view,Qt.Key.Key_Z,Qt.KeyboardModifier.ControlModifier|Qt.KeyboardModifier.ShiftModifier);QTest.qWait(30)
    assert abs(float(mol.frames[0,0,0])-float(coordinates[0,0,0])-.2)<1e-5
    checks.append('real Ctrl+Z and Ctrl+Shift+Z window shortcuts undo and redo without ambiguity')
    QTest.qWait(60);window.grab().save(str(target.with_suffix('.png')))
    report=dict(checks=checks,errors=errors,gpu=view.gpu_info,platform=__import__('sys').platform)
    target.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report));window.modified=False;window.close()


if __name__=='__main__':
    multiprocessing.freeze_support();main()
