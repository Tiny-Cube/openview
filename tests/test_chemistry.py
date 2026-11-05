from pathlib import Path
import numpy as np
import pytest
from openview.formats import load,gaussian_text,save_xyz,BOHR
from openview.model import Molecule,infer_bonds,measure

@pytest.fixture
def water():return Molecule([8,1,1],[[0,0,0],[.9572,0,0],[-.239987,.927297,0]],'Water').infer_bonds()

def write(tmp_path,suffix,text):
    p=tmp_path/('sample'+suffix);p.write_text(text,encoding='utf-8');return p

def test_water_bonds_and_measure(water):
    assert set(map(tuple,water.bonds))=={(0,1,1),(0,2,1)}
    assert '104.51' in measure(water.frames[0],[1,0,2])
    assert '0.9572' in measure(water.frames[0],[0,1])

def test_spatial_hash_matches_distance_reference():
    rng=np.random.default_rng(19);z=rng.choice([1,6,7,8,16,26],300);xyz=rng.uniform(-12,12,(300,3)).astype('f4')
    from openview.model import RADII
    expected={(j,i,1) for i in range(len(z)) for j in range(i) if .16<float(np.sum((xyz[i]-xyz[j])**2))<(1.2*(RADII[z[i]]+RADII[z[j]]))**2}
    assert set(map(tuple,infer_bonds(z,xyz)))==expected

def test_gaussian_roundtrip_and_parity(tmp_path,water):
    result=load(write(tmp_path,'.gjf',gaussian_text(water)))
    np.testing.assert_allclose(result.frames,water.frames,atol=1e-8)
    assert result.charge==0 and result.multiplicity==1
    water.charge=1
    with pytest.raises(ValueError,match='parity'):gaussian_text(water)
    water.multiplicity=2; assert '\n1 2\n' in gaussian_text(water)

def test_gaussian_bohr_input_and_export_route_validation(tmp_path,water):
    text=gaussian_text(water).replace('Opt Freq','Opt Freq Units=Bohr')
    mol=load(write(tmp_path,'.com',text))
    np.testing.assert_allclose(mol.frames,water.frames*BOHR,atol=1e-7)
    assert 'Units=Angstrom' in mol.metadata['route']
    with pytest.raises(ValueError,match='Angstrom'):gaussian_text(water,route='# HF Units=Bohr')
    text=gaussian_text(water,route='# HF Geom=Connectivity')
    assert load(write(tmp_path,'.gjf',text)).bonds.tolist()==[[0,1,1],[0,2,1]]

@pytest.mark.parametrize('route',['# HF Geom=AllCheck','# HF Geom=(Check,Connectivity)','# HF Geom=ZMatrix'])
def test_cartesian_export_rejects_incompatible_geometry(water,route):
    with pytest.raises(ValueError,match='Cartesian export'):
        gaussian_text(water,route=route)

def test_spin_cannot_exceed_electrons():
    hydrogen=Molecule([1],[[0,0,0]],multiplicity=4)
    with pytest.raises(ValueError,match='spin exceeds'):
        gaussian_text(hydrogen)

def test_xyz_trajectory_roundtrip(tmp_path,water):
    water.frames=np.concatenate([water.frames,water.frames+.2])
    p=tmp_path/'中文 molecule.xyz';save_xyz(water,p);result=load(p)
    np.testing.assert_allclose(result.frames,water.frames,atol=1e-7)
    with pytest.raises(ValueError,match='identities'):load(write(tmp_path,'.xyz','1\nA\nH 0 0 0\n1\nB\nHe 1 0 0\n'))

@pytest.mark.parametrize('text',['2\nshort\nH 0 0 0\n','1\nbad\nUnobtainium 0 0 0\n','1\nbad\nH nan 0 0\n'])
def test_invalid_xyz_is_rejected(tmp_path,text):
    with pytest.raises(ValueError):load(write(tmp_path,'.xyz',text))

def orientation(x=0):
    return f''' Standard orientation:
 ---------------------------------------------------------------------
 Center     Atomic      Atomic             Coordinates (Angstroms)
 Number     Number       Type             X           Y           Z
 ---------------------------------------------------------------------
 1 8 0 {x} 0.0 0.0
 2 1 0 0.9572 0.0 0.0
 3 1 0 -0.239987 0.927297 0.0
 ---------------------------------------------------------------------
'''

def test_log_geometry_trajectory_modes_energy(tmp_path):
    text=' Charge = 0 Multiplicity = 1\n'+orientation()+orientation()+orientation(.1)+''' SCF Done: E(RB3LYP) = -76.123456D+00 A.U.
 Frequencies -- -10.0 1000.0 1600.0
 Red. masses -- 1 1 1
 IR Inten -- 2 3 4
 Atom AN X Y Z X Y Z X Y Z
 1 8 .1 0 0 0 .1 0 0 0 .1
 2 1 0 0 0 0 0 0 0 0 0
 3 1 0 0 0 0 0 0 0 0 0
 Normal termination of Gaussian
'''
    mol=load(write(tmp_path,'.log',text))
    assert mol.frames.shape==(2,3,3) and mol.modes.shape==(3,3,3)
    assert mol.frequencies==[-10,1000,1600] and mol.intensities==[2,3,4]
    assert mol.energies==[-76.123456] and mol.metadata['normal_termination']

def test_pdb_models_altloc_and_elements(tmp_path):
    def atom(serial,name,symbol,x,alt=' '):
        return f'HETATM{serial:5d} {name:4s}{alt}MOL A   1    {x:8.3f}{0:8.3f}{0:8.3f}  1.00  0.00          {symbol:>2s}\n'
    text='MODEL        1\n'+atom(1,' CA ','C',0)+atom(2,' CA ','C',.1,'B')+atom(3,' O  ','O',1.3)+'ENDMDL\nMODEL        2\n'+atom(1,' CA ','C',.2)+atom(3,' O  ','O',1.4)+'ENDMDL\n'
    mol=load(write(tmp_path,'.pdb',text));assert mol.numbers.tolist()==[6,8] and mol.frames.shape==(2,2,3)

def test_mol_connectivity_and_charge(tmp_path):
    text='Water\nOpenView\n\n  3  2  0  0  0  0            999 V2000\n'
    for sym,x,y in [('O',0,0),('H',.95,0),('H',-.24,.93)]:text+=f'{x:10.4f}{y:10.4f}{0:10.4f} {sym:<3} 0  0  0  0  0  0  0  0  0  0  0  0\n'
    text+='  1  2  1  0  0  0  0\n  1  3  1  0  0  0  0\nM  END\n$$$$\n'
    mol=load(write(tmp_path,'.sdf',text));assert mol.bonds.tolist()==[[0,1,1],[0,2,1]]

def test_fchk_units(tmp_path):
    text='Hydrogen\nSP RHF STO-3G\n'
    for name,tail in [('Charge','I              0'),('Multiplicity','I              1'),('Atomic numbers','I   N=           2')]:text+=f'{name:<43}{tail}\n'
    text+=' 1 1\n'+f'{"Current cartesian coordinates":<43}R   N=           6\n 0 0 0 1.4 0 0\n'
    mol=load(write(tmp_path,'.fchk',text));assert abs(mol.frames[0,1,0]-1.4*BOHR)<1e-7

@pytest.mark.parametrize('count,unit',[(2,BOHR),(-2,1.)])
def test_cube_grid_and_units(tmp_path,count,unit):
    text=f'Cube\nGenerated test\n1 0 0 0\n{count} 1 0 0\n{count} 0 1 0\n{count} 0 0 1\n1 0 1 0 0\n0 1 2 3 4 5 6 7\n'
    mol=load(write(tmp_path,'.cube',text));assert mol.grid['values'].shape==(2,2,2)
    np.testing.assert_allclose(mol.grid['axes'],np.eye(3)*unit);assert abs(mol.frames[0,0,0]-unit)<1e-7

def test_degenerate_measurement_is_reported():
    assert 'undefined' in measure(np.zeros((4,3)),[0,1,2])
    assert 'undefined' in measure(np.zeros((4,3)),[0,1,2,3])

def test_surface_uses_grid_transform():
    from openview.work import surfaces
    grid=np.indices((12,12,12)).astype('f4')
    values=3.-np.sqrt(np.sum((grid-5.5)**2,axis=0))
    mol=Molecule([1],[[0,0,0]],grid={'values':values,'axes':np.diag([2.,1.,.5]),'origin':np.array([10,0,0])})
    meshes=surfaces(mol,.3);assert len(meshes)==2
    for vertices,normals,color in meshes:
        assert len(vertices)%3==0 and vertices[:,0].min()>10
        np.testing.assert_allclose(np.linalg.norm(normals,axis=1),1,atol=1e-6)
