"""Non-carbon presets, chemical joins, symbol mapping and bounded editing."""
import os
os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
from time import perf_counter
import numpy as np
import pytest
from PySide6.QtWidgets import QApplication
from openview.app import MainWindow
from openview.fragments import CATALOG, FRAGMENTS, available_fragments, build_element, element_options
from openview.formats import gaussian_text, load
from openview.model import Molecule, SYMBOLS


def make(number, key=None):
    return build_element([], np.empty((1, 0, 3)), [], key or element_options(number)['default'],
                         atomic_number=number, position=[0, 0, 0])


@pytest.mark.parametrize('number,hydrogens,key', [(7,3,'tetsss'), (8,2,'bent'), (9,1,'sing'),
    (15,3,'tetsss'), (16,2,'bent'), (17,1,'sing'), (26,6,'oct')])
def test_defaults_and_bond_geometry(number, hydrogens, key):
    assert element_options(number)['default'] == key
    n,p,b,c = make(number)
    assert n.tolist() == [number] + [1]*hydrogens
    assert c == 0 and len(b) == hydrogens
    radius = CATALOG['covalent_radii'][SYMBOLS[number]] + .30
    np.testing.assert_allclose(np.linalg.norm(p[0,1:]-p[0,0], axis=1), radius, atol=1e-6)
    if number == 7:
        assert np.linalg.matrix_rank(p[0,1:]-p[0,0]) == 3  # pyramidal, not planar NH3
    if number == 8:
        assert np.linalg.matrix_rank(p[0,1:]-p[0,0]) == 2  # bent, not linear H2O


def test_all_original_real_element_presets_build_and_modern_symbols_are_not_special_atoms():
    assert len(CATALOG['elements']) == 111 and len(CATALOG['templates']) == 29
    for number in range(1,119):
        options = available_fragments(number)
        assert element_options(number)['default'] in [f.key for f in options]
        if number >= 112:
            assert [f.key for f in options] == ['atom']
            assert SYMBOLS[number] not in CATALOG['elements']
        for option in options:
            n,p,b,c = make(number,option.key)
            assert n[c] == number and np.isfinite(p).all()
            assert len(n) == len(option.orders)+1
            assert sorted(b[:,2]) == sorted(option.orders)
            assert np.count_nonzero(n==0) == sum(order != 1 for order in option.orders)


@pytest.mark.parametrize('number,key,formula,length', [
    (7,'tetsss',{6:1,7:1,1:5},1.47), (8,'bent',{6:1,8:1,1:4},1.43),
    (9,'sing',{6:1,9:1,1:3},1.35)])
def test_heteroatom_attach_and_roundtrip(number,key,formula,length,tmp_path):
    n,p,b,c = make(6)
    n,p,b,c = build_element(n,p,b,key,atomic_number=number,target=1)
    assert dict(zip(*np.unique(n,return_counts=True))) == formula
    assert np.linalg.norm(p[0,c]-p[0,0]) == pytest.approx(length)
    path=tmp_path/'heteroatom.gjf'
    path.write_text(gaussian_text(Molecule(n,p,'Built molecule',bonds=b),include_connectivity=True))
    reopened=load(path)
    np.testing.assert_array_equal(reopened.numbers,n)
    np.testing.assert_allclose(reopened.frames,p,atol=1e-7)
    assert sorted(map(tuple,reopened.bonds))==sorted(map(tuple,b))


def test_nitrogen_triple_attachment_consumes_open_valence():
    n,p,b,c=make(6,'st')
    target=int(b[b[:,2]==3][0,1])
    n,p,b,c=build_element(n,p,b,'dit',atomic_number=7,target=target,match_orders=True)
    assert sorted(n)==[1,6,7]
    assert sorted(b[:,2])==[1,3]
    assert np.linalg.norm(p[0,c]-p[0,0])==pytest.approx(1.47*.78)


def test_ten_coordinate_replacement_is_bounded_and_preserves_heavy_neighbors():
    directions=np.asarray(FRAGMENTS['pentprism'].directions)
    n=np.array([26]+[6]*10)
    p=np.vstack(([0,0,0],directions*2))[None,:,:]
    p=np.concatenate((p,p+[3,4,5]))
    b=np.array([[0,i,1] for i in range(1,11)])
    start=perf_counter()
    n2,p2,b2,c=build_element(n,p,b,'pentprism',atomic_number=27,target=0)
    assert perf_counter()-start < 3  # catches factorial permutation enumeration
    assert n2.tolist()==[27]+[6]*10 and c==0
    np.testing.assert_array_equal(p2,p.astype('f4'))
    np.testing.assert_array_equal(b2,b)


def test_switching_elements_remembers_choices_and_undo_restores_topology():
    app=QApplication.instance() or QApplication([])
    window=MainWindow()
    try:
        n,p,b,_=make(6)
        window.accept_molecule(Molecule(n,p,'Methane',bonds=b))
        original=window.geometry_snapshot()
        window.add_element.setCurrentText('O')
        assert window.element_fragment=='bent'
        assert set(window.fragment_palette.buttons)=={'atom','trid','bent'}
        window.fragment_palette.buttons['bent'].click()
        window.replace_atom(1)
        assert sorted(window.molecule.numbers)==[1,1,1,1,6,8]
        window.undo()
        for before,after in zip(original[:3],window.geometry_snapshot()[:3]):
            np.testing.assert_array_equal(before,after)
        window.redo();assert 8 in window.molecule.numbers
        window.fragment_palette.buttons['trid'].click()
        window.add_element.setCurrentText('N');assert window.element_fragment=='tetsss'
        window.add_element.setCurrentText('O');assert window.element_fragment=='trid'
        window.add_element.setCurrentText('Fe');assert len(window.fragment_palette.buttons)==18
        window.add_element.setCurrentText('Cn');assert list(window.fragment_palette.buttons)==['atom']
        assert 'no original coordination preset' in window.fragment_palette.caption.text()
    finally:
        window.modified=False;window.close()
