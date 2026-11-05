"""Editing data regressions; GPU/pixel checks live in tools/check_interaction.py."""
import os
import time
import numpy as np
import pytest

os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
from PySide6.QtCore import QTimer, Qt
from PySide6.QtTest import QTest
from PySide6.QtWidgets import QApplication, QFileDialog
from openview.app import MainWindow
from openview.formats import load
from openview.model import Molecule


@pytest.fixture(scope='module')
def app():
    return QApplication.instance() or QApplication([])


@pytest.fixture
def window(app):
    window = MainWindow()
    yield window
    for task in list(window.tasks):
        task.cancel(); task.wait(3000)
    window.modified = False
    window.close()


def test_rebuild_bonds_is_undoable_and_invalidates_results(window):
    molecule = Molecule([6, 6], [[0, 0, 0], [1.4, 0, 0]], bonds=[[0, 1, 2]], energies=[-1.])
    window.accept_molecule(molecule)
    window.rebuild_bonds()
    deadline = time.monotonic() + 15
    while window.tasks and time.monotonic() < deadline:
        QTest.qWait(20)
    assert not window.tasks
    assert molecule.bonds[0, 2] == 1 and window.modified
    assert molecule.energies == []
    window.undo()
    assert molecule.bonds[0, 2] == 2
    window.redo()
    assert molecule.bonds[0, 2] == 1


def test_selection_is_not_limited_to_measurement_count(window):
    window.select_all()
    assert len(window.selection) == 12
    assert len(window.table.selectionModel().selectedRows()) == 12
    assert '12 atoms selected' in window.measurement.text()
    window.table.clearSelection(); window.table.selectAll()
    assert len(window.selection) == 12
    window.selection = list(range(5)); window.sync_selection()
    window.delete_atoms()
    assert len(window.molecule.numbers) == 7
    window.undo()
    assert len(window.molecule.numbers) == 12


def test_coordinate_edit_rejects_overflow_and_redo_branch_is_cleared(window):
    before = window.molecule.frames.copy()
    index = window.table_model.index(0, 2)
    assert not window.table_model.setData(index, '1e100')
    np.testing.assert_array_equal(window.molecule.frames, before)
    window.edit_atom(0, 2, 2.)
    window.undo(); window.redo()
    assert window.molecule.frames[0, 0, 0] == 2.
    window.undo(); window.edit_atom(0, 2, 3.)
    assert not window.redo_history


def test_replacing_element_preserves_connections_and_is_undoable(window):
    old = window.molecule.numbers.copy(); bonds = window.molecule.bonds.copy()
    window.add_element.setCurrentText('N'); window.choose_element_fragment('atom'); window.replace_atom(0)
    assert window.molecule.numbers[0] == 7
    np.testing.assert_array_equal(window.molecule.bonds, bonds)
    window.undo(); np.testing.assert_array_equal(window.molecule.numbers, old)


def test_xyz_export_does_not_claim_connectivity_saved(window, tmp_path, monkeypatch):
    window.edit_atom(0, 2, 2.)
    path = tmp_path / 'export.xyz'
    monkeypatch.setattr(QFileDialog, 'getSaveFileName', lambda *a: (str(path), ''))
    window.export_xyz()
    assert path.is_file() and window.modified


def test_save_gjf_clears_dirty_and_roundtrips_bonds(window, app, tmp_path, monkeypatch):
    window.edit_atom(0, 2, 2.)
    path = tmp_path / 'saved.gjf'
    monkeypatch.setattr(QFileDialog, 'getSaveFileName', lambda *a: (str(path), ''))
    QTimer.singleShot(0, lambda: app.activeModalWidget().save_input())
    window.gaussian()
    assert not window.modified
    saved = load(path)
    np.testing.assert_allclose(saved.frames, window.molecule.frames, atol=1e-7)
    assert sorted(map(tuple,saved.bonds)) == sorted(map(tuple,window.molecule.bonds))


def test_undo_stops_playback(window):
    window.edit_atom(0, 2, 2.)
    window.animation.start()
    window.undo()
    assert not window.animation.isActive() and not window.mode_active


def test_alternative_gaussian_charge_does_not_change_working_molecule(window,app,tmp_path,monkeypatch):
    window.edit_atom(0,2,2.)
    path=tmp_path/'cation.gjf'
    monkeypatch.setattr(QFileDialog,'getSaveFileName',lambda *a:(str(path),''))
    def save_cation():
        dialog=app.activeModalWidget();dialog.charge.setValue(1);dialog.spin.setValue(2);dialog.save_input()
    QTimer.singleShot(0,save_cation);window.gaussian()
    assert (window.molecule.charge,window.molecule.multiplicity)==(0,1)
    assert window.modified
    exported=load(path);assert (exported.charge,exported.multiplicity)==(1,2)


@pytest.mark.parametrize('accept',[False,True])
def test_geometry_dialog_preview_cancel_commit(window,app,accept):
    mol=Molecule([6,6],[[0,0,0],[1.4,0,0]],bonds=[[0,1,2]])
    window.accept_molecule(mol);window.selection=[0,1];window.sync_selection()
    before=mol.frames.copy()
    observed=[]
    def operate():
        dialog=app.activeModalWidget();dialog.value.setValue(1.8)
        observed.append(window.view.positions[1,0])
        np.testing.assert_array_equal(mol.frames,before)
        dialog.accept() if accept else dialog.reject()
    QTimer.singleShot(0,operate);window.edit_geometry()
    assert observed[0]==pytest.approx(1.8)
    assert mol.frames[0,1,0]==pytest.approx(1.8 if accept else 1.4)
    assert window.modified==accept and len(window.history)==int(accept)
    assert mol.bonds[0,2]==2
    if accept:
        window.undo();np.testing.assert_array_equal(mol.frames,before)
