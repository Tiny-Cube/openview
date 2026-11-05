"""One owned worker process per task; Qt never parses or waits on the pipe."""
import multiprocessing as mp
import traceback
import numpy as np
from PySide6.QtCore import QThread, Signal

def surfaces(mol,level):
    from skimage.measure import marching_cubes
    grid=mol.grid; values=grid['values']; meshes=[]
    axes=grid['axes']
    if abs(np.linalg.det(axes))<1e-15: raise ValueError('CUBE grid axes are singular')
    for sign,color in ((1,(.3,.83,.66)),(-1,(.65,.48,.96))):
        value=sign*level
        if not float(values.min())<value<float(values.max()): continue
        vertices,faces,normals,_=marching_cubes(values,value,allow_degenerate=False)
        world=vertices@axes+grid['origin']
        world_normals=normals@np.linalg.inv(axes).T
        world_normals/=np.maximum(np.linalg.norm(world_normals,axis=1,keepdims=True),1e-15)
        meshes.append((world[faces].reshape(-1,3).astype('f4'),world_normals[faces].reshape(-1,3).astype('f4'),color))
    if not meshes: raise ValueError('Isovalue lies outside the scalar grid range')
    return meshes

def worker(conn,kind,args):
    try:
        if kind=='load':
            from .formats import load
            result=load(*args)
        elif kind=='bonds':
            from .model import infer_bonds
            result=infer_bonds(*args)
        elif kind=='surface': result=surfaces(*args)
        else: raise ValueError('Unknown worker operation')
        conn.send((True,result))
    except Exception as exc: conn.send((False,f'{type(exc).__name__}: {exc}'))
    finally: conn.close()

class Task(QThread):
    completed=Signal(object)
    failed=Signal(str)

    def __init__(self,kind,*args,parent=None):
        super().__init__(parent); self.kind=kind; self.args=args; self.process=None

    def run(self):
        ctx=mp.get_context('spawn'); recv,send=ctx.Pipe(duplex=False)
        process=ctx.Process(target=worker,args=(send,self.kind,self.args),daemon=True)
        self.process=process
        try:
            process.start(); send.close()
            while not self.isInterruptionRequested():
                if recv.poll(.05):
                    success,result=recv.recv()
                    if not self.isInterruptionRequested():
                        if success: self.completed.emit(result)
                        else: self.failed.emit(result)
                    break
                if not process.is_alive():
                    if not self.isInterruptionRequested(): self.failed.emit(f'Worker exited with code {process.exitcode}')
                    break
        except Exception as exc:
            if not self.isInterruptionRequested(): self.failed.emit(str(exc))
        finally:
            recv.close(); send.close()
            if process.pid is not None:
                if process.is_alive(): process.terminate()
                process.join(timeout=2)
                if process.is_alive(): process.kill(); process.join()
                process.close()
            self.process=None

    def cancel(self): self.requestInterruption()
