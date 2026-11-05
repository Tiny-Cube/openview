"""Chemical topology, geometry, serialization and undo for carbon building."""
import os
os.environ.setdefault('QT_QPA_PLATFORM','offscreen')
import numpy as np
import pytest
from PySide6.QtWidgets import QApplication
from openview.fragments import build_carbon, FRAGMENTS
from openview.model import Molecule
from openview.formats import gaussian_text, load, save_xyz
from openview.app import MainWindow


def fragment(key):
    n,p,b,_=build_carbon([],np.empty((1,0,3)),[],key,position=[0,0,0])
    return Molecule(n,p,key,bonds=b)


@pytest.mark.parametrize('key,atoms,ports,orders',[
    ('atom',1,0,[]),('st',3,1,[1,3]),('dd',3,2,[2,2]),
    ('ssd',4,1,[1,1,2]),('saa',4,2,[1,1.5,1.5]),('tet',5,0,[1,1,1,1])])
def test_coordination_templates(key,atoms,ports,orders):
    mol=fragment(key)
    assert len(mol.numbers)==atoms and np.count_nonzero(mol.numbers==0)==ports
    assert sorted(mol.bonds[:,2])==orders
    vectors=mol.frames[0,1:]-mol.frames[0,0]
    if len(vectors):
        vectors/=np.linalg.norm(vectors,axis=1)[:,None]
        dots=(vectors@vectors.T)[np.triu_indices(len(vectors),1)]
        expected=-1 if key in ('st','dd') else -1/3 if key=='tet' else -.5
        np.testing.assert_allclose(dots,expected,atol=1e-6)


@pytest.mark.parametrize('key,order,formula,length',[
    ('tet',1,{6:2,1:6},1.54),('ssd',2,{6:2,1:4},1.3552),('st',3,{6:2,1:2},1.2012)])
def test_build_ethane_ethene_ethyne_and_roundtrip(key,order,formula,length,tmp_path):
    mol=fragment(key);target=int(mol.bonds[mol.bonds[:,2]==order][0,1])
    before=[a.copy() for a in (mol.numbers,mol.frames,mol.bonds)]
    n,p,b,center=build_carbon(mol.numbers,mol.frames,mol.bonds,key,target=target)
    for old,actual in zip(before,(mol.numbers,mol.frames,mol.bonds)):np.testing.assert_array_equal(old,actual)
    assert dict(zip(*np.unique(n,return_counts=True)))==formula
    carbon=np.flatnonzero(n==6)
    assert np.linalg.norm(p[0,carbon[1]]-p[0,carbon[0]])==pytest.approx(length)
    assert [bond[2] for bond in b if set(bond[:2])==set(carbon)]==[order]
    if key=='ssd':assert np.linalg.matrix_rank(p[0]-p[0,0],tol=1e-5)==2
    complete=Molecule(n,p,'Built molecule',bonds=b)
    path=tmp_path/'built.gjf';path.write_text(gaussian_text(complete,include_connectivity=True))
    reopened=load(path)
    np.testing.assert_array_equal(reopened.numbers,n);np.testing.assert_allclose(reopened.frames,p,atol=1e-7)
    assert sorted(map(tuple,reopened.bonds))==sorted(map(tuple,b))


def test_mismatch_is_transactional_and_open_valences_not_exported(tmp_path):
    mol=fragment('st');before=mol.frames.copy()
    with pytest.raises(ValueError,match='matching free bond'):
        build_carbon(mol.numbers,mol.frames,mol.bonds,'tet',target=2,match_orders=True)
    np.testing.assert_array_equal(mol.frames,before)
    with pytest.raises(ValueError,match='open-valence'):gaussian_text(mol)
    path=tmp_path/'unfinished.xyz'
    with pytest.raises(ValueError,match='open-valence'):save_xyz(mol,path)
    assert not path.exists()


def test_replacement_retains_heavy_bonds_and_each_trajectory_frame():
    mol=fragment('tet')
    n,p,b,c=build_carbon(mol.numbers,mol.frames,mol.bonds,'tet',target=1)
    p=np.concatenate((p,p+np.array([[[3,4,5]]],dtype='f4')))
    n2,p2,b2,c2=build_carbon(n,p,b,'ssd',target=0)
    assert len(n2)==7 and np.count_nonzero(n2==0)==1
    np.testing.assert_allclose(p2[1]-p2[0],np.tile([3,4,5],(len(n2),1)),atol=1e-6)
    np.testing.assert_array_equal(p2[:,c2],p[:,0])
    assert any(np.array_equal(p2[:,i],p[:,c]) for i in np.flatnonzero(n2==6))
    assert sorted(b2[:,2])==[1,1,1,1,1,2]


def test_palette_actions_undo_and_failed_build_does_not_clear_history():
    app=QApplication.instance() or QApplication([]);window=MainWindow()
    try:
        window.accept_molecule(fragment('tet'))
        window.carbon_palette.buttons['tet'].click()
        assert window.view.interaction_mode=='add'
        original=window.geometry_snapshot()
        window.replace_atom(1)
        assert np.count_nonzero(window.molecule.numbers==6)==2 and len(window.history)==1
        window.undo()
        for expected,actual in zip(original[:3],window.geometry_snapshot()[:3]):np.testing.assert_array_equal(actual,expected)
        window.redo();assert len(window.molecule.numbers)==8
        window.accept_molecule(Molecule([6,6,6],[[0,0,0],[1,0,0],[-1,0,0]],bonds=[[0,1,1],[0,2,1]]))
        # Zero available ports is the deliberate bare-atom mode, so use a
        # three-heavy-neighbor target with the two-coordinate linear pattern.
        window.accept_molecule(Molecule([6,6,6,6],[[0,0,0],[1,0,0],[-1,0,0],[0,1,0]],bonds=[[0,1,1],[0,2,1],[0,3,1]]))
        window.choose_carbon_fragment('dd');window.replace_atom(0)
        assert not window.history and len(window.molecule.numbers)==4
        assert 'matching free bond' in window.status.text()
        window.choose_carbon_fragment('dd');window.place_atom([5,0,0])
        assert window.molecule.numbers[-3:].tolist()==[6,0,0]
        assert window.table_model.data(window.table_model.index(len(window.molecule.numbers)-1,1))=='Open'
    finally:
        window.modified=False;window.close()


def test_default_join_preserves_target_order_even_if_template_differs():
    mol=fragment('st')
    n,p,b,center=build_carbon(mol.numbers,mol.frames,mol.bonds,'tet',target=2)
    assert [order for a,c,order in b if {int(a),int(c)}=={0,center}]==[3]


def test_alkene_planarity_does_not_depend_on_camera_orientation():
    from openview.fragments import axial_rotation
    mol=fragment('ssd')
    for angle in (.3,1.4,2.9):
        n,p,b,c=build_carbon(mol.numbers,mol.frames,mol.bonds,'ssd',target=3,
                           orientation=axial_rotation([.2,.5,.8],angle))
        assert np.linalg.matrix_rank(p[0]-p[0,0],tol=1e-5)==2
