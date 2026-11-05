"""Readers for documented chemistry formats. Never execute content from input files."""
from pathlib import Path
import re
import numpy as np
from .model import Molecule, SYMBOLS, element
from .checkpoint import BOHR, bond_table, read_chk
from .zmatrix import build_zmatrix

FLOAT = r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[EeDd][-+]?\d+)?"

def real(s): return float(s.replace('D','E').replace('d','e'))

def read_xyz(path):
    frames=[]; numbers=None; title=''
    with open(path,encoding='utf-8-sig',errors='replace') as f:
        while True:
            line=f.readline()
            if not line: break
            if not line.strip(): continue
            try: n=int(line)
            except ValueError: raise ValueError('XYZ frame must start with an atom count')
            if n<1 or n>2_000_000: raise ValueError('XYZ atom count out of range')
            comment=f.readline().strip(); zs=[]; xyz=[]
            for _ in range(n):
                row=f.readline().split()
                if len(row)<4: raise ValueError('Truncated XYZ coordinates')
                zs.append(element(row[0])); xyz.append([real(x) for x in row[1:4]])
            if numbers is None: numbers=zs; title=comment
            elif zs!=numbers: raise ValueError('XYZ trajectory changes atom identities or count')
            frames.append(xyz)
    return Molecule(numbers or [],frames,title or Path(path).stem)

def read_gaussian_input(path):
    lines=Path(path).read_text(encoding='utf-8-sig',errors='replace').splitlines()
    link1=next((i for i,line in enumerate(lines) if line.strip().lower()=='--link1--'),len(lines))
    lines=lines[:link1]
    start=next((i for i,s in enumerate(lines) if s.lstrip().startswith('#')),None)
    if start is None: raise ValueError('Missing Gaussian route section (#)')
    i=start; route=[]
    while i<len(lines) and lines[i].strip(): route.append(lines[i].strip()); i+=1
    while i<len(lines) and not lines[i].strip(): i+=1
    title=[]
    while i<len(lines) and lines[i].strip(): title.append(lines[i]); i+=1
    while i<len(lines) and not lines[i].strip(): i+=1
    if i>=len(lines): raise ValueError('Missing charge and multiplicity')
    cm=lines[i].split()
    if len(cm)!=2: raise ValueError('Fragment-specific charge/multiplicity inputs are not yet supported')
    charge,mult=int(cm[0]),int(cm[1]); i+=1
    numbers=[]; xyz=[];rows=[]
    while i<len(lines) and lines[i].strip():
        row=lines[i].replace(',',' ').split(); i+=1
        rows.append(row)
    if not rows:raise ValueError('No explicit molecule in GJF; checkpoint-only geometry needs its referenced CHK file')
    is_zmatrix=len(rows[0])==1
    if is_zmatrix:numbers,xyz=build_zmatrix(rows,lines[i:])
    else:
        for row in rows:
            if len(row) not in (4,5):raise ValueError('Invalid Cartesian Gaussian atom row')
            symbol=row[0].split('(')[0].split('-')[0]
            numbers.append(element(symbol)); xyz.append([real(v) for v in row[-3:]])
    # Explicit connectivity following Cartesian geometry, if supplied.
    bonds=[]
    for s in (lines[i:] if re.search(r'connectivity',' '.join(route),re.I) else []):
        row=s.split()
        if len(row)<3 or not row[0].isdigit(): continue
        try:
            a=int(row[0])-1
            for j in range(1,len(row)-1,2): bonds.append((a,int(row[j])-1,float(row[j+1])))
        except ValueError: continue
    route_text=' '.join(route)
    if re.search(r'\bunits\s*=\s*\(?\s*(?:bohr|au|a0)\b',route_text,re.I):
        xyz=(np.asarray(xyz)*BOHR).tolist()
        route_text=re.sub(r'\bunits\s*=\s*\(?\s*(?:bohr|au|a0)\s*\)?','Units=Angstrom',route_text,flags=re.I)
    return Molecule(numbers,xyz,' '.join(title),bonds=bonds if re.search(r'connectivity',route_text,re.I) else None,charge=charge,multiplicity=mult,metadata={'route':route_text,'input_coordinates':'z-matrix' if is_zmatrix else 'cartesian'})

def read_gaussian_log(path):
    frames=[]; numbers=None; energies=[]; frequencies=[]; intensities=[]; modes=[]
    charge=0; mult=1; normal=False
    # Avoid alternating input/standard orientations for the same optimization step.
    with open(path,encoding='utf-8',errors='replace') as probe:
        orientation='Standard orientation:' if any('Standard orientation:' in s for s in probe) else 'Input orientation:'
    with open(path,encoding='utf-8',errors='replace') as f:
        iterator=iter(f)
        for line in iterator:
            if orientation in line:
                # The five-line Gaussian table heading ends at the second rule.
                rules=0
                for s in iterator:
                    if s.strip().startswith('-----'): rules+=1
                    if rules==2: break
                zs=[]; xyz=[]
                for s in iterator:
                    if s.strip().startswith('-----'): break
                    row=s.split()
                    if len(row)>=6:
                        zs.append(int(row[1])); xyz.append([real(v) for v in row[-3:]])
                if zs:
                    if numbers is None: numbers=zs
                    if zs!=numbers: raise ValueError('Multi-job log changes molecule; split jobs before opening')
                    coords=np.asarray(xyz,dtype='f4')
                    if not frames or not np.allclose(frames[-1],coords,atol=1e-6,rtol=0): frames.append(coords)
            elif 'SCF Done:' in line:
                m=re.search(r'=\s*('+FLOAT+')',line)
                if m: energies.append(real(m[1]))
            elif 'Charge =' in line and 'Multiplicity =' in line:
                m=re.search(r'Charge\s*=\s*(-?\d+)\s+Multiplicity\s*=\s*(\d+)',line)
                if m: charge,mult=map(int,m.groups())
            elif 'Frequencies --' in line:
                fs=[real(x) for x in line.split('--',1)[1].split()]; frequencies.extend(fs)
                # Follow this mode group through the displacement table.
                found=False
                for s in iterator:
                    if 'IR Inten' in s: intensities.extend(real(x) for x in s.split('--',1)[1].split())
                    if re.match(r'\s*Atom\s+AN\s+',s): found=True; break
                    if len(s)>0 and 'Normal termination' in s: normal=True; break
                if found and numbers:
                    block=[]
                    for _ in numbers:
                        row=next(iterator,'').split()
                        if len(row)!=2+3*len(fs): raise ValueError('Incomplete vibration displacement table')
                        block.append([real(x) for x in row[2:]])
                    data=np.asarray(block,dtype='f4').reshape(len(numbers),len(fs),3)
                    modes.extend(data[:,k,:] for k in range(len(fs)))
            if 'Normal termination of Gaussian' in line: normal=True
    if not frames: raise ValueError('No Cartesian orientation tables found in Gaussian log')
    return Molecule(numbers,frames,Path(path).stem,energies=energies,frequencies=frequencies,
                    intensities=intensities,modes=np.asarray(modes) if modes else None,
                    charge=charge,multiplicity=mult,metadata={'normal_termination':normal})

def read_pdb(path):
    frames=[]; numbers=[]; bonds=set(); serials={}; coords=[]; first=True
    with open(path,encoding='utf-8',errors='replace') as f:
        for line in f:
            record=line[:6].strip()
            if record in ('ATOM','HETATM'):
                if line[16:17] not in ('',' ','A'): continue
                name=line[12:16]; symbol=line[76:78].strip()
                if not symbol:
                    clean=re.sub(r'[^A-Za-z]','',name)
                    # PDB right-aligned atom names distinguish carbon CA from calcium CA.
                    symbol=clean[:1] if name.startswith(' ') else clean[:2]
                    if symbol.lower() not in {s.lower() for s in SYMBOLS}: symbol=clean[:1]
                if first:
                    serials[int(line[6:11])]=len(numbers); numbers.append(element(symbol))
                elif len(coords)>=len(numbers) or element(symbol)!=numbers[len(coords)]: raise ValueError('PDB MODEL atom identity mismatch')
                coords.append([real(line[30:38]),real(line[38:46]),real(line[46:54])])
            elif record=='ENDMDL':
                if coords: frames.append(coords); coords=[]; first=False
            elif record=='CONECT':
                ids=[int(line[i:i+5]) for i in range(6,len(line.rstrip()),5) if line[i:i+5].strip()]
                for other in ids[1:]:
                    if ids[0] in serials and other in serials and ids[0]!=other: bonds.add(tuple(sorted((serials[ids[0]],serials[other])))+(1,))
    if coords: frames.append(coords)
    # PDB connectivity can be incomplete (ligands only); infer missing pairs too.
    mol=Molecule(numbers,frames,Path(path).stem)
    mol.infer_bonds()
    if bonds: mol.bonds=np.asarray(sorted(set(map(tuple,mol.bonds))|bonds),dtype='i4')
    return mol

def read_mol(path):
    lines=Path(path).read_text(encoding='utf-8-sig',errors='replace').splitlines()
    if len(lines)<4 or 'V3000' in lines[3]: raise ValueError('Expected MDL V2000 MOL/SDF (first record)')
    n,b=int(lines[3][:3]),int(lines[3][3:6]); zs=[]; xyz=[]; bonds=[]
    for line in lines[4:4+n]:
        xyz.append([real(line[i:i+10]) for i in (0,10,20)]); zs.append(element(line[31:34]))
    if len(xyz)!=n: raise ValueError('Truncated MOL atoms')
    for line in lines[4+n:4+n+b]: bonds.append((int(line[:3])-1,int(line[3:6])-1,int(line[6:9])))
    charge=0
    for line in lines[4+n+b:]:
        if line.startswith('M  CHG'): charge+=sum(map(int,line.split()[4::2]))
    return Molecule(zs,xyz,lines[0],bonds=bonds,charge=charge)

def read_fchk(path):
    fields={}
    selected={'Atomic numbers','Current cartesian coordinates','Charge','Multiplicity','Total Energy',
              'Number of atoms','Number of electrons','Number of basis functions','Full Title','Route',
              'Number of Normal Modes','Vib-E2','Vib-Modes','MxBond','NBond','IBond','RBond','Mulliken Charges'}
    header=re.compile(r'^(.{43})([IRCLH])\s+(?:(N=)\s*(\d+)|(.*?))\s*$')
    with open(path,encoding='utf-8',errors='replace') as f:
        title=f.readline().strip(); next(f,'')
        for line in f:
            match=header.match(line.rstrip('\r\n'))
            if not match:
                if line.strip():raise ValueError('Invalid FCHK field header: '+line[:80].strip())
                continue
            raw_label,kind,array,count,scalar=match.groups();label=raw_label.strip();wanted=label in selected
            if array:
                count=int(count);values=[];consumed=0
                if kind=='C':
                    rows=[]
                    for _ in range((count+5)//6):
                        row=f.readline()
                        if not row:raise ValueError(f'Truncated FCHK character array: {label}')
                        rows.append(row.rstrip('\r\n'))
                    if wanted:fields[label]=''.join(rows)[:count*12].strip()
                    continue
                while consumed<count:
                    row=f.readline()
                    if not row or header.match(row.rstrip('\r\n')):raise ValueError(f'Truncated FCHK array: {label}')
                    tokens=list(row.strip()) if kind in ('H','L') else re.findall(FLOAT,row)
                    if not tokens or consumed+len(tokens)>count:raise ValueError(f'Invalid FCHK array length: {label}')
                    consumed+=len(tokens)
                    if wanted:values.extend(real(v) for v in tokens)
                if wanted:fields[label]=values
            elif wanted:fields[label]=scalar.strip() if kind=='C' else real(scalar.strip())
    if 'Atomic numbers' not in fields or 'Current cartesian coordinates' not in fields: raise ValueError('FCHK has no Cartesian geometry')
    n=len(fields['Atomic numbers'])
    if 'Number of atoms' in fields and int(fields['Number of atoms'])!=n:raise ValueError('FCHK atom count mismatch')
    xyz=np.asarray(fields['Current cartesian coordinates']).reshape(-1,3)*BOHR
    bonds=None;metadata={'checkpoint_reader':'formatted'}
    if {'MxBond','NBond','IBond','RBond'}<=fields.keys():
        bonds=bond_table(n,fields['MxBond'],fields['NBond'],fields['IBond'],fields['RBond']);metadata['bond_source']='FCHK bond table'
    for key,label in [('route','Route'),('number_of_electrons','Number of electrons'),('basis_functions','Number of basis functions'),('mulliken_charges','Mulliken Charges')]:
        if label in fields:metadata[key]=fields[label]
    frequencies=[];intensities=[];modes=None;nm=int(fields.get('Number of Normal Modes',0))
    if nm and 'Vib-E2' in fields:
        vib=np.asarray(fields['Vib-E2'])
        if len(vib)<4*nm:raise ValueError('Incomplete FCHK vibration properties')
        frequencies=vib[:nm].tolist();intensities=vib[3*nm:4*nm].tolist()
    if nm and 'Vib-Modes' in fields:
        if len(fields['Vib-Modes'])!=nm*n*3:raise ValueError('FCHK vibration displacement dimensions mismatch')
        modes=np.asarray(fields['Vib-Modes'],dtype='f4').reshape(nm,n,3)
    return Molecule(fields['Atomic numbers'],xyz,fields.get('Full Title',title),bonds=bonds,
                    charge=int(fields.get('Charge',0)),multiplicity=int(fields.get('Multiplicity',1)),
                    energies=[fields['Total Energy']] if 'Total Energy' in fields else [],
                    frequencies=frequencies,intensities=intensities,modes=modes,metadata=metadata)

def read_cube(path):
    with open(path,encoding='utf-8',errors='replace') as f:
        title=f.readline().strip(); next(f)
        row=f.readline().split(); n=int(row[0]); origin=np.array([real(v) for v in row[1:4]])
        if n<0: raise ValueError('Multi-orbital CUBE is not supported; export one orbital per file')
        shape=[]; axes=[]; signs=[]
        for _ in range(3):
            row=f.readline().split(); count=int(row[0]); signs.append(np.sign(count)); shape.append(abs(count)); axes.append([real(x) for x in row[1:4]])
        if len(set(signs))!=1 or 0 in signs: raise ValueError('Inconsistent CUBE axis units/counts')
        unit=BOHR if signs[0]>0 else 1.
        zs=[]; xyz=[]
        for _ in range(n):
            row=f.readline().split(); zs.append(element(row[0])); xyz.append([real(x)*unit for x in row[2:5]])
        values=np.fromstring(f.read().replace('D','E'),sep=' ',dtype='f4')
        if len(values)!=int(np.prod(shape)): raise ValueError('CUBE scalar grid size mismatch')
    return Molecule(zs,xyz,title,grid={'origin':origin*unit,'axes':np.asarray(axes)*unit,'values':values.reshape(shape)})

READERS={'.xyz':read_xyz,'.gjf':read_gaussian_input,'.com':read_gaussian_input,'.log':read_gaussian_log,'.out':read_gaussian_log,'.pdb':read_pdb,'.ent':read_pdb,'.mol':read_mol,'.sdf':read_mol,'.fchk':read_fchk,'.fch':read_fchk,'.chk':read_chk,'.cube':read_cube,'.cub':read_cube}

def load(path):
    reader=READERS.get(Path(path).suffix.lower())
    if not reader: raise ValueError('Supported formats: '+', '.join(READERS))
    mol=reader(path).infer_bonds(); mol.metadata['source']=str(path)
    return mol

def gaussian_text(mol,frame=0,route='#p B3LYP/6-31G(d) Opt Freq',memory='2GB',cores=4,checkpoint='molecule.chk',include_connectivity=False):
    if np.any(mol.numbers==0): raise ValueError('Connect or delete the orange open-valence handles before exporting Gaussian input.')
    if mol.multiplicity<1: raise ValueError('Multiplicity must be positive')
    electrons=int(mol.numbers.sum())-mol.charge
    if electrons<0 or mol.multiplicity-1>electrons or (electrons-(mol.multiplicity-1))%2: raise ValueError('Charge and spin multiplicity have incompatible electron parity or spin exceeds electron count')
    if any('\n' in str(v) or '\r' in str(v) for v in (memory,checkpoint)): raise ValueError('Link0 values must be single-line')
    route=route.strip()
    if re.search(r'\n\s*\n',route): raise ValueError('Route section must not contain a blank line')
    geometry=re.search(r'\bgeom\s*=\s*(\([^)]*\)|[^\s]+)',route,re.I)
    if geometry and re.search(r'\b(?:allcheck|check|checkpoint|zmatrix)\b',geometry[1],re.I):
        raise ValueError('Remove checkpoint/Z-matrix Geom options before Cartesian export')
    # Coordinates emitted below are always Cartesian Angstroms.
    if re.search(r'\bunits\s*=\s*\(?\s*(?:bohr|au|a0)\b',route,re.I): raise ValueError('Export coordinates use Angstroms; change Units=Bohr to Units=Angstrom')
    include_connectivity=include_connectivity or bool(re.search(r'\bconnectivity\b',route,re.I))
    if include_connectivity and not re.search(r'\bconnectivity\b',route,re.I):
        geom=re.search(r'\bgeom\s*=\s*(\([^)]*\)|[^\s]+)',route,re.I)
        if geom:
            options=geom[1].strip('()')
            if re.search(r'\b(?:allcheck|check|zmatrix)\b',options,re.I):raise ValueError('Remove checkpoint/Z-matrix Geom options before Cartesian export')
            route=route[:geom.start()]+'Geom=('+options+',Connectivity)'+route[geom.end():]
        else:route+=' Geom=Connectivity'
    if not route.startswith('#'): route='#p '+route
    title=mol.title.replace('\n',' ').replace('\r',' ') or 'OpenView molecule'
    head=f'%chk={checkpoint}\n%mem={memory}\n%nprocshared={int(cores)}\n{route}\n\n{title}\n\n{mol.charge} {mol.multiplicity}\n'
    output=head+'\n'.join(f'{SYMBOLS[z]:<3} {p[0]: .9f} {p[1]: .9f} {p[2]: .9f}' for z,p in zip(mol.numbers,mol.frames[frame]))+'\n\n'
    if include_connectivity:
        neighbors={i:[] for i in range(len(mol.numbers))}
        for a,b,order in mol.bonds if mol.bonds is not None else []:
            a,b=sorted((int(a),int(b)));neighbors[a].append((b,float(order)))
        output+='\n'.join(str(i+1)+''.join(f' {b+1} {order:g}' for b,order in sorted(neighbors[i])) for i in neighbors)+'\n\n'
    return output

def save_xyz(mol,path,frame=None):
    if np.any(mol.numbers==0): raise ValueError('Connect or delete the orange open-valence handles before exporting XYZ.')
    frames=mol.frames if frame is None else mol.frames[frame:frame+1]
    with open(path,'w',encoding='utf-8',newline='\n') as f:
        for xyz in frames:
            f.write(f'{len(mol.numbers)}\n{mol.title}\n')
            for z,p in zip(mol.numbers,xyz): f.write(f'{SYMBOLS[z]} {p[0]:.8f} {p[1]:.8f} {p[2]:.8f}\n')
