import struct
import numpy as np
import pytest
from openview.formats import load,gaussian_text
from openview.checkpoint import read_native_chk,UnsupportedCheckpoint,BOHR
from openview.model import measure

def binary_chk(endian='<'):
    """Independent tiny CHK fixture with relocated, reordered logical records."""
    n=3;stride=n+3;core=np.zeros(20+5*stride,dtype=endian+'i8')
    core[:10]=[-12349,20,n,0,1,5,5,10,7,7];core[20:20+n]=[8,1,1];core[20+n:20+stride]=-3
    xyz=np.array([[0,0,0],[1.8,0,0],[-.45,1.75,0]],dtype=endian+'f8')
    floating=core.view(endian+'f8');floating[20+stride:20+stride+n]=[8,1,1];floating[20+2*stride:20+2*stride+3*n]=xyz.ravel()
    bonds=struct.pack(endian+'2i',n,2)+np.array([2,1,1,2,3,1,0,1,0],dtype=endian+'i8').tobytes()+np.array([1.5,1.,1.5,0,1,0],dtype=endian+'f8').tobytes()
    blocks=[(698,np.concatenate([xyz.ravel(),xyz.ravel()]).astype(endian+'f8').tobytes()),(997,core.tobytes()),(503,bonds)]
    data_start=64;offset=data_start;rows=[];content=bytearray(data_start*8)
    for identifier,block in blocks:
        words=len(block)//8;rows.extend([identifier,offset,offset+words,offset+words]);content.extend(block);offset+=words
    content[:32]=struct.pack(endian+'4q',len(blocks),offset,data_start,0)
    content[32:32+len(rows)*8]=np.array(rows,dtype=endian+'i8').tobytes()
    return content,xyz

@pytest.mark.parametrize('endian',['<','>'])
def test_native_checkpoint_directory_offsets_and_orders(tmp_path,endian):
    data,xyz=binary_chk(endian);p=tmp_path/'standalone.chk';p.write_bytes(data)
    mol=load(p);assert mol.numbers.tolist()==[8,1,1] and mol.charge==0 and mol.multiplicity==1
    np.testing.assert_allclose(mol.frames[0],xyz*BOHR,atol=1e-7)
    assert mol.bonds.tolist()==[[0,1,1.5],[0,2,1.]] and mol.metadata['checkpoint_reader']=='native-64'

def test_bad_directory_is_rejected_and_auxiliary_geometry_is_not_substituted(tmp_path):
    data,_=binary_chk();p=tmp_path/'broken.chk'
    struct.pack_into('<q',data,32+16,len(data)//8+1);p.write_bytes(data)
    with pytest.raises(ValueError,match='out-of-file'):read_native_chk(p)
    data,xyz=binary_chk();struct.pack_into('<d',data,64*8,10.);p.write_bytes(data)
    mol=read_native_chk(p)
    np.testing.assert_allclose(mol.frames[0],xyz*BOHR,atol=1e-7)
    assert mol.metadata['geometry_crosscheck'].startswith('different auxiliary geometry')

def test_native_checkpoint_accepts_standard_orientation(tmp_path):
    data,xyz=binary_chk()
    rotation=np.array([[0.,-1.,0.],[1.,0.,0.],[0.,0.,1.]])
    standard=xyz@rotation+np.array([4.,-2.,1.])
    data[64*8:64*8+xyz.nbytes]=standard.astype('<f8').tobytes()
    p=tmp_path/'rotated.chk';p.write_bytes(data)
    mol=load(p)
    np.testing.assert_allclose(mol.frames[0],xyz*BOHR,atol=1e-7)
    assert mol.metadata['geometry_crosscheck']=='proper rigid alignment'

def test_geometry_crosscheck_rejects_reflection_and_distortion():
    from openview.checkpoint import same_geometry
    xyz=np.array([[0.,0.,0.],[1.,0.,0.],[0.,2.,0.],[0.,0.,3.]])
    assert same_geometry(xyz,xyz+10)
    assert not same_geometry(xyz,xyz*np.array([-1,1,1]))
    distorted=xyz.copy();distorted[0,0]=.01
    assert not same_geometry(xyz,distorted)

def test_unknown_chk_never_uses_same_named_fchk(tmp_path,monkeypatch):
    monkeypatch.delenv('OPENVIEW_FORMCHK',raising=False);monkeypatch.delenv('GVNEXT_FORMCHK',raising=False);monkeypatch.setattr('openview.checkpoint.shutil.which',lambda _:None)
    p=tmp_path/'test.chk';p.write_bytes(b'not a checkpoint')
    p.with_suffix('.fchk').write_text('This must not be substituted')
    with pytest.raises(ValueError,match='No companion file was substituted'):load(p)

def test_truncated_fchk_array_has_clear_failure(tmp_path):
    p=tmp_path/'bad.fchk';p.write_text('Bad\nSP HF STO-3G\n'+f'{"Atomic numbers":<43}I N= 2\n 1\n'+f'{"Charge":<43}I 0\n')
    with pytest.raises(ValueError,match='Truncated FCHK array'):load(p)

def test_zmatrix_variables_and_dihedral(tmp_path):
    p=tmp_path/'internal.gjf';p.write_text('# HF/STO-3G\n\nInternal geometry\n\n0 1\nC\nC 1 R\nH 2 H 1 A\nH 2 H 1 A 3 D\n\nR=1.5\nH=1.1\nA=109.5\nD=60\n')
    mol=load(p);assert mol.numbers.tolist()==[6,6,1,1]
    assert '1.5000' in measure(mol.frames[0],[0,1])
    assert '109.500' in measure(mol.frames[0],[2,1,0])
    assert '60.000' in measure(mol.frames[0],[3,1,0,2])

def test_zmatrix_dummy_centers_and_undefined_variable(tmp_path):
    p=tmp_path/'dummy.gjf';p.write_text('# HF\n\nDummy\n\n0 1\nX\nH 1 1\nH 1 1 2 90\n\n')
    mol=load(p);assert len(mol.numbers)==2
    p.write_text('# HF\n\nUndefined\n\n0 1\nH\nH 1 Missing\n\n')
    with pytest.raises(ValueError,match='Undefined Z-matrix'):load(p)

def test_fractional_connectivity_roundtrips(tmp_path):
    data,_=binary_chk();p=tmp_path/'molecule.chk';p.write_bytes(data);mol=load(p)
    gjf=tmp_path/'edited.gjf';gjf.write_text(gaussian_text(mol,include_connectivity=True))
    roundtrip=load(gjf);np.testing.assert_array_equal(roundtrip.bonds,mol.bonds)
