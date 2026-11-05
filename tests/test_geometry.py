import numpy as np
import pytest
from openview.geometry import geometry_value, edited_geometry


def test_distance_translates_group_without_changing_internal_geometry():
    xyz=np.array([[0.,0.,0.],[1.,0.,0.],[1.,1.,0.]],dtype='f4')
    bonds=np.array([[0,1,2],[1,2,1]])
    result=edited_geometry(xyz,bonds,[0,1],1.5)
    np.testing.assert_allclose(result,[[0,0,0],[1.5,0,0],[1.5,1,0]])
    np.testing.assert_array_equal(xyz,[[0,0,0],[1,0,0],[1,1,0]])
    assert geometry_value(result,[1,2])==1


@pytest.mark.parametrize('target',[30.,109.5,150.])
def test_angle_rotation_keeps_bond_lengths(target):
    xyz=np.array([[1.,0.,0.],[0.,0.,0.],[0.,1.,0.],[0.,2.,1.]])
    bonds=np.array([[0,1,1],[1,2,1],[2,3,1]])
    result=edited_geometry(xyz,bonds,[0,1,2],target)
    assert geometry_value(result,[0,1,2])==pytest.approx(target)
    assert geometry_value(result,[1,2])==pytest.approx(1.)
    assert geometry_value(result,[2,3])==pytest.approx(np.sqrt(2))
    np.testing.assert_array_equal(result[:2],xyz[:2])


@pytest.mark.parametrize('target',[-170.,-60.,0.,60.,170.])
def test_dihedral_rotation_and_fixed_side(target):
    xyz=np.array([[0.,1.,0.],[0.,0.,0.],[1.,0.,0.],[1.,0.,1.],[2.,0.,1.]])
    bonds=np.array([[0,1,1],[1,2,1],[2,3,1],[3,4,1]])
    result=edited_geometry(xyz,bonds,[0,1,2,3],target)
    assert geometry_value(result,[0,1,2,3])==pytest.approx(target,abs=1e-7)
    np.testing.assert_allclose(result[:3],xyz[:3])
    assert geometry_value(result,[3,4])==pytest.approx(1.)


def test_ring_requires_explicit_atom_only_edit():
    xyz=np.array([[0.,0.,0.],[1.,0.,0.],[0.,1.,0.]])
    bonds=np.array([[0,1,1],[1,2,1],[2,0,1]])
    with pytest.raises(ValueError,match='ring'):
        edited_geometry(xyz,bonds,[0,1],1.5)
    result=edited_geometry(xyz,bonds,[0,1],1.5,move_group=False)
    np.testing.assert_array_equal(result[[0,2]],xyz[[0,2]])
    assert geometry_value(result,[0,1])==1.5


def test_undefined_dihedral_is_rejected():
    with pytest.raises(ValueError,match='collinear'):
        geometry_value(np.array([[0.,0.,0.],[1.,0.,0.],[2.,0.,0.],[3.,0.,0.]]),[0,1,2,3])
