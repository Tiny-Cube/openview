"""Exercise non-carbon fragment selection and attachment with real Qt clicks."""
import argparse
import json
import multiprocessing
from pathlib import Path
import numpy as np
from PySide6.QtCore import QPoint, Qt
from PySide6.QtTest import QTest
from PySide6.QtWidgets import QApplication, QToolBar
from openview.app import MainWindow, configure_graphics, STYLE
from openview.formats import gaussian_text
from openview.model import SYMBOLS


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--show',action='store_true')
    args=parser.parse_args();args.out.mkdir(parents=True,exist_ok=True)
    configure_graphics();app=QApplication([]);app.setStyle('Fusion');app.setStyleSheet(STYLE)
    window=MainWindow();window.setWindowTitle('OpenView · Element Fragment verification');window.resize(1380,860)
    window.show();window.activateWindow();app.setActiveWindow(window);QTest.qWait(250);view=window.view
    assert view.ctx is not None and not view.error,view.error
    toolbar=window.findChild(QToolBar)
    new=next(a for a in toolbar.actions() if a.text()=='New molecule')
    checks=[]

    def choose(symbol):
        combo=window.add_element;combo.setFocus()
        QTest.keyClick(combo,Qt.Key.Key_Home)
        for _ in range(combo.findText(symbol)):QTest.keyClick(combo,Qt.Key.Key_Down)
        assert combo.currentText()==symbol
        QTest.qWait(30)

    def fresh():
        window.modified=False
        QTest.mouseClick(toolbar.widgetForAction(new),Qt.MouseButton.LeftButton);QTest.qWait(80)

    def click_atom(index):
        xy,_=view.project()
        QTest.mouseClick(view,Qt.MouseButton.LeftButton,pos=QPoint(*np.rint(xy[index]).astype(int)))
        QTest.qWait(60)

    for symbol,key,size in [('N','tetsss',4),('O','bent',3),('F','sing',2),('P','tetsss',4),('S','bent',3),('Fe','oct',7)]:
        choose(symbol)
        assert window.element_fragment==key
        QTest.mouseClick(window.fragment_palette.buttons[key],Qt.MouseButton.LeftButton)
        fresh()
        assert window.molecule.numbers[0]==SYMBOLS.index(symbol) and len(window.molecule.numbers)==size
        assert not view.grabFramebuffer().isNull()
        if symbol in ('N','O','Fe'):assert window.grab().save(str(args.out/(symbol+'-fragments.png')))
        checks.append(symbol+': keyboard element selection, default palette click, New molecule and OpenGL render')
    bar=window.fragment_palette.scroll.horizontalScrollBar()
    assert bar.maximum()>0
    bar.setValue(bar.maximum());QTest.qWait(50)
    last=window.fragment_palette.buttons['pentprism']
    assert window.fragment_palette.scroll.viewport().rect().contains(last.mapTo(window.fragment_palette.scroll.viewport(),last.rect().center()))
    QTest.mouseClick(last,Qt.MouseButton.LeftButton);fresh()
    assert len(window.molecule.numbers)==11
    checks.append('Fe: horizontal scrolling exposes ten-coordinate template and builds 11 atoms')
    for symbol,key,formula,name in [('N','tetsss',[1,1,1,1,1,6,7],'methylamine'),('O','bent',[1,1,1,1,6,8],'methanol')]:
        choose('C');QTest.mouseClick(window.fragment_palette.buttons['tet'],Qt.MouseButton.LeftButton);fresh()
        choose(symbol);QTest.mouseClick(window.fragment_palette.buttons[key],Qt.MouseButton.LeftButton)
        click_atom(1)
        assert sorted(window.molecule.numbers)==formula
        QTest.keyClick(view,Qt.Key.Key_Z,Qt.KeyboardModifier.ControlModifier);QTest.qWait(30)
        assert sorted(window.molecule.numbers)==[1,1,1,1,6]
        window.redo();assert sorted(window.molecule.numbers)==formula
        (args.out/(name+'.gjf')).write_text(gaussian_text(window.molecule,include_connectivity=True))
        view.fit();QTest.qWait(80)
        assert window.grab().save(str(args.out/(name+'.png')))
        checks.append(name+': viewport attachment consumes H, Ctrl+Z/redo restores topology, GJF export succeeds')
    report={'checks':checks,'gpu':view.gpu_info,'errors':[]}
    (args.out/'element-palette-check.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report),flush=True);window.modified=False
    if args.show:return app.exec()
    window.close()


if __name__=='__main__':multiprocessing.freeze_support();main()
