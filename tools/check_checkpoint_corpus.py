"""Read-only real CHK regression with independent companion LOG comparison."""
import argparse
import hashlib
import json
from pathlib import Path
import time
import numpy as np
from openview.formats import load, gaussian_text
from openview.checkpoint import same_geometry
import tempfile


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('directory',type=Path)
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args();records=[]
    paths=sorted(args.directory.glob('*.chk'))
    if not paths:raise SystemExit('No CHK files found')
    for path in paths:
        started=time.perf_counter();digest=hashlib.sha256(path.read_bytes()).hexdigest()
        record={'file':path.name,'sha256':digest}
        try:
            mol=load(path)
            record.update(atoms=len(mol.numbers),bonds=len(mol.bonds),reader=mol.metadata.get('checkpoint_reader'),
                geometry_source=mol.metadata.get('geometry_source'),crosscheck=mol.metadata.get('geometry_crosscheck'))
            log=path.with_suffix('.log')
            if log.exists():
                reference=load(log)
                np.testing.assert_array_equal(mol.numbers,reference.numbers)
                if not same_geometry(mol.frames[0].astype('f8'),reference.frames[-1].astype('f8'),tolerance=2e-5):
                    raise AssertionError('Molecular geometry disagrees with final companion LOG orientation')
                record['last_log_geometry']='matches under proper rigid transform (2e-5 Å tolerance)'
                record['log_sha256']=hashlib.sha256(log.read_bytes()).hexdigest()
            with tempfile.TemporaryDirectory() as directory:
                exported=Path(directory)/'roundtrip.gjf';exported.write_text(gaussian_text(mol,include_connectivity=True))
                back=load(exported)
                np.testing.assert_allclose(back.frames,mol.frames,atol=1e-6,rtol=0)
                assert sorted(map(tuple,back.bonds))==sorted(map(tuple,mol.bonds))
            assert hashlib.sha256(path.read_bytes()).hexdigest()==digest
            record['status']='passed'
        except Exception as error:
            record.update(status='failed',error=str(error))
        record['seconds']=time.perf_counter()-started;records.append(record)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(records,indent=2))
    print(json.dumps(records,indent=2))
    return int(any(r['status']!='passed' for r in records))


if __name__=='__main__':raise SystemExit(main())
