import argparse
import json
from pathlib import Path
import time
import tempfile
import numpy as np
from openview.formats import load,save_xyz
from openview.model import Molecule

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--out',default='analysis/core-benchmark.json');args=parser.parse_args();results=[]
    for n in (1000,10000,100000):
        # Repeated separated water molecules, with known connectivity count.
        molecules=(n+2)//3;side=int(np.ceil(molecules**(1/3)))
        centers=np.indices((side,side,side)).reshape(3,-1).T[:molecules]*4.
        xyz=(centers[:,None,:]+np.array([[0,0,0],[.9572,0,0],[-.239987,.927297,0]])).reshape(-1,3)[:n].astype('f4')
        z=np.resize([8,1,1],n);mol=Molecule(z,xyz,f'{n} atom benchmark')
        with tempfile.TemporaryDirectory() as folder:
            path=Path(folder)/'molecule.xyz';save_xyz(mol,path);start=time.perf_counter();cpu=time.process_time();parsed=load(path);elapsed=time.perf_counter()-start;cpu=time.process_time()-cpu
        expected=(n//3)*2+max(0,n%3-1)
        assert len(parsed.bonds)==expected,(n,len(parsed.bonds),expected)
        results.append({'atoms':n,'bonds':len(parsed.bonds),'read_and_bond_seconds':elapsed,'cpu_seconds':cpu,'geometry_bytes':parsed.frames.nbytes+parsed.bonds.nbytes+parsed.numbers.nbytes})
    path=Path(args.out);path.parent.mkdir(parents=True,exist_ok=True);path.write_text(json.dumps(results,indent=2),encoding='utf-8');print(json.dumps(results))

if __name__=='__main__':main()
