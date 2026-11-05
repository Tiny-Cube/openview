"""Gaussian checkpoint geometry readers.

Native CHK profile: 64-bit word directory, record 997 magic -12349/header 20.
Offsets are relative to validated logical records, never absolute file offsets.
Unknown profiles fail explicitly or use the operator's Gaussian formchk utility.
"""
from pathlib import Path
import os
import shutil
import struct
import subprocess
import tempfile
import numpy as np
from .model import Molecule

BOHR=.529177210903

class UnsupportedCheckpoint(ValueError):pass

def same_geometry(first, second, tolerance=1e-7):
    """Compare atom-matched geometries allowing a proper rigid transform.

    Gaussian can store input-oriented coordinates in 997 and standard-oriented
    coordinates in 698/701. Reflections and actual geometry changes must fail.
    This uses O(N) storage rather than an all-pairs distance matrix.
    """
    if first.shape != second.shape or not np.isfinite(first).all() or not np.isfinite(second).all():
        return False
    a = first - first.mean(axis=0)
    b = second - second.mean(axis=0)
    u, _, vt = np.linalg.svd(a.T @ b)
    correction = np.eye(3)
    correction[-1, -1] = np.linalg.det(u @ vt)
    return bool(np.max(np.abs(a @ (u @ correction @ vt) - b)) <= tolerance)

def bond_table(n,maximum,counts,indices,orders):
    maximum=int(maximum);counts=np.asarray(counts,dtype='i8')
    indices=np.asarray(indices,dtype='i8').reshape(n,maximum)
    orders=np.asarray(orders,dtype='f8').reshape(n,maximum)
    if counts.shape!=(n,) or np.any(counts<0) or np.any(counts>maximum):raise ValueError('Invalid checkpoint bond counts')
    pairs={}
    for a in range(n):
        for j in range(int(counts[a])):
            b=int(indices[a,j])-1;order=float(orders[a,j])
            if b<0 or b>=n or b==a or not np.isfinite(order) or order<=0:raise ValueError('Invalid checkpoint bond table')
            key=tuple(sorted((a,b)))
            if key in pairs and abs(pairs[key]-order)>1e-6:raise ValueError('Inconsistent reciprocal checkpoint bond orders')
            pairs[key]=order
    return np.array([(*pair,order) for pair,order in sorted(pairs.items())],dtype='f4').reshape(-1,3)

def read_native_chk(path):
    path=Path(path);size=path.stat().st_size
    with path.open('rb') as f:
        header=f.read(32);profile=None
        if len(header)!=32:raise UnsupportedCheckpoint('Truncated CHK directory header')
        for endian in ('<','>'):
            count,total,data_start,reserved=struct.unpack(endian+'4q',header)
            if 0<count and (4+4*count)*8<=data_start*8<=size and total*8==size and reserved==0:
                profile=(endian,count,data_start);break
        if profile is None:raise UnsupportedCheckpoint('Unknown CHK word size, byte order or directory layout')
        endian,count,data_start=profile
        raw=f.read(count*32)
        rows=np.frombuffer(raw,dtype=endian+'i8').reshape(count,4)
        records={};intervals=[]
        for identifier,start,end,allocated in rows:
            identifier,start,end,allocated=map(int,(identifier,start,end,allocated))
            if not data_start<=start<=end<=allocated<=size//8:raise ValueError('CHK directory contains an out-of-file record')
            if identifier:
                if identifier<0 or identifier in records:raise ValueError('Invalid/duplicate CHK logical record')
                records[identifier]=(start,end);intervals.append((start,allocated))
        intervals.sort()
        if any(b>c for (_,b),(c,_) in zip(intervals,intervals[1:])):raise ValueError('Overlapping CHK logical records')
        def record(identifier):
            if identifier not in records:raise UnsupportedCheckpoint(f'CHK logical record {identifier} is absent')
            start,end=records[identifier];f.seek(start*8);return f.read((end-start)*8)
        data=record(997);ints=np.frombuffer(data,dtype=endian+'i8');floats=np.frombuffer(data,dtype=endian+'f8')
        if len(ints)<20 or tuple(ints[:2])!=(-12349,20):raise UnsupportedCheckpoint('Unrecognized CHK molecular record version')
        n,charge,mult=map(int,ints[2:5]);stride=n+3;coordinate_start=20+2*stride
        if n<1 or coordinate_start+3*n>len(ints) or mult<1:raise ValueError('Invalid CHK molecular dimensions')
        z=ints[20:20+n].astype('i4')
        if np.any(z<1) or np.any(z>118) or not np.all(ints[20+n:20+stride]==-3):raise UnsupportedCheckpoint('Unrecognized CHK center table layout')
        if int(ints[5]+ints[6])!=int(ints[7]) or int(ints[5]-ints[6])!=mult-1:raise ValueError('Inconsistent CHK electron/spin fields')
        nuclear=floats[20+stride:20+stride+n]
        if not np.isfinite(nuclear).all() or np.any(nuclear<0) or np.any(nuclear>z+1e-6):raise ValueError('Invalid CHK nuclear charges')
        xyz=floats[coordinate_start:coordinate_start+3*n].reshape(n,3).copy()
        metadata={'checkpoint_reader':'native-64','byte_order':'little' if endian=='<' else 'big','logical_records':sorted(records),
                  'number_of_electrons':int(ints[7]),'basis_functions':int(ints[8]),'native_profile':'997/-12349/20'}
        # 997 is the molecular geometry. 698 can contain a standard-oriented
        # geometry from an earlier step (observed in real excited/anion jobs).
        # Its disagreement alone does not establish corruption of record 997.
        if 698 in records:
            check=np.frombuffer(record(698),dtype=endian+'f8')
            if len(check)==6*n:
                metadata['geometry_crosscheck_record']=698
                metadata['geometry_crosscheck']=('proper rigid alignment' if same_geometry(check[:3*n].reshape(n,3),xyz)
                    else 'different auxiliary geometry; molecular record 997 retained')
        metadata['geometry_source']='CHK molecular record 997'
        title=path.stem
        if 502 in records:
            text=record(502)
            if len(text)>=112 and struct.unpack_from(endian+'q',text)[0]==-12345:
                candidate=text[32:112].decode('ascii',errors='replace').strip()
                if candidate and '\ufffd' not in candidate:title=candidate
        bonds=None
        if 503 in records:
            block=record(503)
            if len(block)>=8:
                bn,maximum=struct.unpack_from(endian+'2i',block)
                if bn==n and maximum>0 and len(block)==8*(1+n+2*n*maximum):
                    bi=np.frombuffer(block,dtype=endian+'i8');br=np.frombuffer(block,dtype=endian+'f8')
                    bonds=bond_table(n,maximum,bi[1:1+n],bi[1+n:1+n+n*maximum],br[1+n+n*maximum:])
                    metadata['bond_source']='CHK logical record 503'
        energies=[]
        if 701 in records:
            block=record(701);bi=np.frombuffer(block,dtype=endian+'i8');br=np.frombuffer(block,dtype=endian+'f8')
            if len(bi)>=13+3*n and bi[4]==3*n and same_geometry(br[13:13+3*n].reshape(n,3),xyz) and np.isfinite(br[9]):
                energies=[float(br[9])];metadata['energy_source']='CHK logical record 701'
        return Molecule(z,xyz*BOHR,title,bonds=bonds,charge=charge,multiplicity=mult,energies=energies,metadata=metadata)

def read_chk(path):
    try:return read_native_chk(path)
    except UnsupportedCheckpoint as native_error:
        configured=os.environ.get('OPENVIEW_FORMCHK') or os.environ.get('GVNEXT_FORMCHK')
        converter=configured or shutil.which('formchk')
        if not converter:
            raise ValueError(f'{native_error}. This CHK variant needs Gaussian formchk: set OPENVIEW_FORMCHK to its executable, or convert to .fchk. No companion file was substituted.') from native_error
        from .formats import read_fchk
        with tempfile.TemporaryDirectory(prefix='openview-formchk-') as folder:
            output=Path(folder)/'converted.fchk'
            with (Path(folder)/'diagnostics.txt').open('w+b') as diagnostics:
                try:
                    result=subprocess.run([str(converter),str(Path(path).resolve()),str(output)],stdin=subprocess.DEVNULL,stdout=diagnostics,stderr=subprocess.STDOUT,timeout=90,check=False)
                except (OSError,subprocess.TimeoutExpired) as exc:raise ValueError(f'formchk conversion failed: {exc}') from exc
                diagnostics.seek(0);message=diagnostics.read(4096).decode(errors='replace')
                if result.returncode or not output.is_file():raise ValueError(f'formchk failed (exit {result.returncode}): {message}')
            mol=read_fchk(output);mol.metadata['checkpoint_reader']='formchk';return mol
