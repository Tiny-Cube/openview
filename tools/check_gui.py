"""Real Qt/OpenGL integration and frame benchmarks; run under xvfb-run on Linux."""
import argparse
import json
import multiprocessing
import os
from pathlib import Path
import sys
import time
import numpy as np
from PySide6.QtCore import QTimer,Qt
from PySide6.QtWidgets import QApplication,QMessageBox
from openview.app import MainWindow,configure_graphics,STYLE
from openview.model import Molecule

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--atoms',type=int,default=10000);parser.add_argument('--frames',type=int,default=180);parser.add_argument('--out',default='analysis/gui-check.json');parser.add_argument('--trajectory',action='store_true');parser.add_argument('--full-bonds',action='store_true');args=parser.parse_args()
    configure_graphics();app=QApplication([]);app.setStyle('Fusion');app.setStyleSheet(STYLE)
    window=MainWindow();window.resize(1280,800);window.show();errors=[];result={};ticks=[];rendered=[];started=time.perf_counter();state={'frame':0,'stage':'init','idle':0}
    window.view.force_full_bonds=args.full_bonds;result['full_detail_bonds']=args.full_bonds
    # Convert unexpected GUI warnings to evidence instead of blocking CI forever.
    QMessageBox.warning=lambda parent,title,text:errors.append(title+': '+text)
    beat=QTimer();beat.setTimerType(Qt.TimerType.PreciseTimer);beat.setInterval(10);beat.timeout.connect(lambda:ticks.append(time.perf_counter()));beat.start()
    spin=QTimer();spin.setTimerType(Qt.TimerType.PreciseTimer);spin.setInterval(0)
    def fail(message):
        errors.append(message);finish()
    def finish():
        spin.stop();beat.stop();window.animation.stop()
        for task in window.tasks:task.cancel()
        for task in window.tasks:task.wait(3000)
        window.modified=False
        result['errors']=errors;result['gpu']=window.view.gpu_info;result['atoms']=args.atoms;result['trajectory']=args.trajectory;result['platform']=sys.platform
        result['viewport_pixels']=[round(window.view.width()*window.devicePixelRatioF()),round(window.view.height()*window.devicePixelRatioF())]
        if len(rendered)>30:
            intervals=np.diff(rendered[20:])*1000
            result['completed_frame_interval_ms']={'median':float(np.median(intervals)),'p95':float(np.percentile(intervals,95)),'max':float(max(intervals))};result['fps']=1000/float(np.mean(intervals))
        if len(ticks)>2:
            intervals=np.diff(ticks)*1000;result['ui_heartbeat_ms']={'median':float(np.median(intervals)),'p95':float(np.percentile(intervals,95)),'max':float(max(intervals))}
        times=window.view.paint_times[20:]
        if times:result['cpu_paint_ms']={'median':float(np.median(times)),'p95':float(np.percentile(times,95))}
        result['draw_calls']=window.view.draw_calls;result['elapsed_seconds']=time.perf_counter()-started
        path=Path(args.out);path.parent.mkdir(parents=True,exist_ok=True)
        # Save a frame with the actual molecule, before the worker replaces it.
        window.grab().save(str(path.with_suffix('.png')))
        path.write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(result));app.exit(1 if errors else 0)
    def rotated():
        if state['stage']!='benchmark':return
        rendered.append(time.perf_counter());state['frame']+=1
        if state['frame']>=args.frames:
            state['stage']='idle';spin.stop();before=len(window.view.paint_times)
            def idle_check():
                result['idle_repaint_count']=len(window.view.paint_times)-before
                window.grab().save(str(Path(args.out).with_suffix('.benchmark.png')))
                # A real worker process must open a file while Qt heartbeats continue.
                state['stage']='worker';state['worker_ticks']=len(ticks);state['worker_started']=time.perf_counter()
                window.open_path(Path(__file__).resolve().parents[1]/'examples'/'water.xyz')
                wait.start()
            QTimer.singleShot(400,idle_check)
    def frame():
        if state['stage']!='benchmark':return
        window.view.yaw+=.012
        if args.trajectory:
            mol=window.molecule
            window.view.set_positions(mol.frames[0]+np.array([0,0,np.sin(state['frame']*.08)*.2],dtype='f4'))
        window.view.update()
    spin.timeout.connect(frame);window.view.frameSwapped.connect(rotated)
    def begin():
        if not window.view.ctx:fail(window.view.error or 'No GL context');return
        side=int(np.ceil(args.atoms**(1/3)));coords=np.indices((side,side,side)).reshape(3,-1).T[:args.atoms].astype('f4')*1.45
        z=np.resize(np.array([6,6,6,8,7,1],dtype='i4'),args.atoms)
        # Known chain connectivity isolates rendering from the parser benchmark.
        bonds=np.column_stack((np.arange(args.atoms-1),np.arange(1,args.atoms),np.ones(args.atoms-1,dtype='i4')))
        mol=Molecule(z,coords,f'GPU benchmark · {args.atoms:,} atoms',bonds=bonds)
        window.accept_molecule(mol);window.view.paint_times.clear();state['stage']='benchmark';spin.start()
    wait=QTimer();wait.setInterval(30)
    def loaded():
        if window.loading:return
        wait.stop()
        result['background_load_seconds']=time.perf_counter()-state['worker_started'];result['heartbeats_during_load']=len(ticks)-state['worker_ticks']
        if len(window.molecule.numbers)!=3:errors.append('Worker did not return water')
        window.pick(0);window.pick(1)
        if '0.9572' not in window.measurement.text():errors.append('Selection measurement mismatch')
        window.edit_atom(1,2,1.1)
        state['stage']='edit';state['deadline']=time.perf_counter()+10
        edit_wait.start()
    wait.timeout.connect(loaded)
    edit_wait=QTimer();edit_wait.setInterval(30)
    def edited():
        if window.tasks and time.perf_counter()<state['deadline']:return
        edit_wait.stop();window.undo()
        if abs(float(window.molecule.frames[0,1,0])-.9572)>1e-5:errors.append('Coordinate undo failed')
        window.selection=[1];window.delete_atoms();window.undo()
        if len(window.molecule.numbers)!=3:errors.append('Topology undo failed')
        window.modified=False
        result['integration_checks']=['GPU rendering','background file load','atom selection','distance','coordinate edit/undo','delete/undo']
        finish()
    edit_wait.timeout.connect(edited)
    QTimer.singleShot(250,begin);QTimer.singleShot(90000,lambda:fail('Timed out') if state['stage']!='done' else None)
    return app.exec()

if __name__=='__main__':
    multiprocessing.freeze_support();raise SystemExit(main())
