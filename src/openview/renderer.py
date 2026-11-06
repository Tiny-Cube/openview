"""Qt6 OpenGL 3.3 instanced renderer. Idle views do not repaint."""
import math
import time
import numpy as np
import moderngl
from PySide6.QtCore import Qt, Signal, QRect
from PySide6.QtGui import QPainter, QColor
from PySide6.QtOpenGLWidgets import QOpenGLWidget
from PySide6.QtWidgets import QRubberBand
from .model import COLORS, RADII

SPHERE_VERT='''#version 330
in vec2 corner;
in vec3 position;
in float radius;
in vec3 color;
uniform mat3 rotation;
uniform vec3 origin;
uniform vec3 extent;
out vec2 uv;
out vec3 tint;
out float center_z;
out float rad;
void main(){
 vec3 p=rotation*(position-origin);
 uv=corner; tint=color; center_z=p.z; rad=radius;
 gl_Position=vec4((p.xy+corner*radius)/extent.xy,-p.z/extent.z,1.0);
}'''
SPHERE_FRAG='''#version 330
in vec2 uv;
in vec3 tint;
in float center_z;
in float rad;
uniform vec3 extent;
out vec4 frag;
void main(){
 float q=dot(uv,uv); if(q>1.0) discard;
 vec3 n=vec3(uv,sqrt(1.0-q));
 float light=.28+.72*max(dot(n,normalize(vec3(-.4,.65,1))),0.0);
 float spec=pow(max(dot(n,normalize(vec3(-.2,.3,1))),0.0),42.0)*.3;
 frag=vec4(tint*light+spec,1);
 gl_FragDepth=.5-(center_z+n.z*rad)/(2.0*extent.z);
}'''
BOND_VERT='''#version 330
in vec3 vertex;
in vec3 start;
in vec3 finish;
in vec3 color_a;
in vec3 color_b;
uniform mat3 rotation;
uniform vec3 origin;
uniform vec3 extent;
uniform float bond_radius;
out vec3 normal;
flat out vec3 tint_a;
flat out vec3 tint_b;
out float along;
void main(){
 vec3 a=rotation*(start-origin),b=rotation*(finish-origin);
 vec3 d=normalize(b-a);
 vec3 u=normalize(cross(d,abs(d.z)<.9?vec3(0,0,1):vec3(0,1,0)));
 vec3 v=cross(d,u);
 normal=vertex.x*u+vertex.y*v;
 vec3 p=mix(a,b,vertex.z)+bond_radius*normal;
 tint_a=color_a; tint_b=color_b; along=vertex.z;
 gl_Position=vec4(p.xy/extent.xy,-p.z/extent.z,1);
}'''
SURFACE_VERT='''#version 330
in vec3 position;
in vec3 normal_in;
uniform mat3 rotation;
uniform vec3 origin;
uniform vec3 extent;
out vec3 normal;
void main(){ vec3 p=rotation*(position-origin); normal=rotation*normal_in;
gl_Position=vec4(p.xy/extent.xy,-p.z/extent.z,1); }'''
LIT_FRAG='''#version 330
in vec3 normal;
flat in vec3 tint_a;
flat in vec3 tint_b;
in float along;
out vec4 frag;
void main(){vec3 tint=along<.5?tint_a:tint_b;float l=.3+.7*abs(dot(normalize(normal),normalize(vec3(-.4,.65,1))));frag=vec4(tint*l,1);}'''
SURFACE_FRAG='''#version 330
in vec3 normal;
uniform vec3 tint;
out vec4 frag;
void main(){float l=.25+.75*abs(dot(normalize(normal),normalize(vec3(-.4,.65,1))));frag=vec4(tint*l,1);}'''

class MolecularView(QOpenGLWidget):
    atomPicked=Signal(int)
    selectionRequested=Signal(object,bool)
    placementRequested=Signal(object)
    replacementRequested=Signal(int)
    selectAllRequested=Signal()
    clearSelectionRequested=Signal()
    deleteRequested=Signal()
    exitPlacement=Signal()
    stats=Signal(str)
    graphicsError=Signal(str)

    def __init__(self,parent=None):
        super().__init__(parent)
        self.setMinimumSize(280,100); self.setFocusPolicy(Qt.FocusPolicy.StrongFocus)
        self.ctx=None; self.error=''; self.molecule=None; self.positions=None
        self.origin=np.zeros(3,dtype='f4'); self.yaw=.3; self.pitch=-.25; self.roll=0.; self.scale=5.; self.depth=20.
        self.pan=np.zeros(2,dtype='f4'); self.style='Ball & stick'; self.selection=[]
        self.dirty=True; self.bonds_dirty=True; self.surface_dirty=True; self.surfaces=[]; self.surface_resources=[]
        self.atom_resources=[]; self.bond_resources=[]; self.last=None; self.dragged=False
        self.paint_times=[]; self.gpu_info=''; self.draw_calls=0
        self.interaction_mode='select'
        self.click_select_held=False
        self.fragment_select_held=False
        self.selection_mode='click'; self.held_selection_modes=[]
        self.selection_drag=None; self.brush_seen=set()
        self.rubber_band=QRubberBand(QRubberBand.Shape.Rectangle,self)
        self.rubber_band.setAttribute(Qt.WidgetAttribute.WA_TransparentForMouseEvents)
        self.force_full_bonds=False

    def initializeGL(self):
        try:
            self.ctx=moderngl.create_context(require=330)
            self.sphere=self.ctx.program(vertex_shader=SPHERE_VERT,fragment_shader=SPHERE_FRAG)
            self.bond=self.ctx.program(vertex_shader=BOND_VERT,fragment_shader=LIT_FRAG)
            self.surface=self.ctx.program(vertex_shader=SURFACE_VERT,fragment_shader=SURFACE_FRAG)
            self.quad=self.ctx.buffer(np.array([[-1,-1],[1,-1],[-1,1],[1,1]],dtype='f4').tobytes())
            def cylinder_buffer(sides):
                vertices=[]
                for j in range(sides):
                    a,b=2*math.pi*j/sides,2*math.pi*(j+1)/sides
                    for low,high in ((0,.5),(.5,1)):
                        vertices.extend([(math.cos(a),math.sin(a),low),(math.cos(b),math.sin(b),low),(math.cos(a),math.sin(a),high),
                                         (math.cos(a),math.sin(a),high),(math.cos(b),math.sin(b),low),(math.cos(b),math.sin(b),high)])
                return self.ctx.buffer(np.asarray(vertices,dtype='f4').tobytes())
            self.cylinder=cylinder_buffer(10)
            self.cylinder_coarse=cylinder_buffer(4)
            self.gpu_info=self.ctx.info['GL_RENDERER']+' / '+self.ctx.info['GL_VERSION']
            self.context().aboutToBeDestroyed.connect(self.cleanup)
            self.dirty=True;self.bonds_dirty=True
        except Exception as exc:
            self.error=f'OpenGL 3.3 initialization failed: {exc}'
            self.graphicsError.emit(self.error)

    def cleanup(self):
        if not self.ctx: return
        self.makeCurrent()
        for resource in self.atom_resources+self.bond_resources+self.surface_resources:
            resource.release()
        for name in ('quad','cylinder','cylinder_coarse','sphere','bond','surface'):
            obj=getattr(self,name,None)
            if obj: obj.release()
        self.atom_resources=[]; self.bond_resources=[]; self.surface_resources=[]
        self.ctx=None; self.doneCurrent()

    def rotation(self):
        cy,sy=math.cos(self.yaw),math.sin(self.yaw); cp,sp=math.cos(self.pitch),math.sin(self.pitch)
        cr,sr=math.cos(self.roll),math.sin(self.roll)
        roll=np.array([[cr,-sr,0],[sr,cr,0],[0,0,1]],dtype='f4')
        return roll@np.array([[cy,0,sy],[sp*sy,cp,-sp*cy],[-cp*sy,sp,cp*cy]],dtype='f4')

    def extent(self):
        aspect=max(1,self.width())/max(1,self.height())
        return np.array([self.scale*aspect,self.scale,max(self.scale*20,self.depth)],dtype='f4')

    def set_molecule(self,mol):
        self.cancel_selection_drag()
        self.molecule=mol; self.positions=mol.frames[0]; self.selection=[]
        self.surfaces=[]; self.surface_dirty=True; self.dirty=True;self.bonds_dirty=True; self.fit()

    def fit(self):
        if self.positions is None: return
        self.origin=((self.positions.min(axis=0)+self.positions.max(axis=0))*.5).astype('f4')
        radius=float(np.max(np.linalg.norm(self.positions-self.origin,axis=1)))+1.5
        self.scale=max(2.,radius*1.12/min(1.,max(1,self.width())/max(1,self.height())))
        self.depth=max(20.,radius*3); self.pan[:]=0; self.update()

    def set_positions(self,positions):
        self.positions=positions; self.dirty=True;self.bonds_dirty=True; self.update()

    def set_style(self,style): self.style=style; self.dirty=True;self.bonds_dirty=True; self.update()
    def set_selection(self,selection): self.selection=list(selection); self.dirty=True; self.update()

    def set_interaction_mode(self,mode):
        self.cancel_selection_drag()
        self.interaction_mode=mode
        self.setCursor(Qt.CursorShape.CrossCursor if mode=='add' else Qt.CursorShape.ArrowCursor)

    def set_selection_mode(self,mode):
        if mode not in ('click','rectangle','brush'):
            raise ValueError('Unknown selection mode')
        self.cancel_selection_drag(); self.selection_mode=mode

    def cancel_selection_drag(self):
        self.selection_drag=None; self.brush_seen.clear(); self.last=None
        self.rubber_band.hide()

    def effective_selection_mode(self):
        if self.click_select_held or self.fragment_select_held:return 'click'
        if self.held_selection_modes:return self.held_selection_modes[-1][1]
        return self.selection_mode if self.interaction_mode=='select' else 'click'

    def hit_atom(self,position):
        if self.positions is None or not len(self.positions):return -1
        xy,z=self.project(); point=np.array([position.x(),position.y()])
        distance=np.linalg.norm(xy-point,axis=1)
        radii=np.clip(RADII[self.molecule.numbers]*.48,.19,.8)
        if self.style=='Space filling':radii=np.clip(RADII[self.molecule.numbers]*1.65,.55,2.8)
        if self.style=='Sticks':radii=np.full(len(z),.14)
        screen=radii*self.height()/(2*self.scale)
        hit=np.flatnonzero(distance<np.maximum(7,screen))
        if not len(hit):return -1
        distance_world=distance[hit]*2*self.scale/max(1,self.height())
        front=z[hit]+np.sqrt(np.maximum(0,radii[hit]**2-distance_world**2))
        return int(hit[np.argmax(front)])

    def brush_at(self,position):
        index=self.hit_atom(position)
        if index>=0 and index not in self.brush_seen:
            self.brush_seen.add(index); self.selectionRequested.emit([index],True)

    def _upload(self):
        if self.molecule is None: return
        z=self.molecule.numbers
        radii=np.clip(RADII[z]*.48,.19,.8)
        if self.style=='Space filling': radii=np.clip(RADII[z]*1.65,.55,2.8)
        if self.style=='Sticks': radii=np.full(len(z),.14)
        # Selection is an atom-only visual property, never a bond endpoint color.
        element_colors=COLORS[z]
        colors=element_colors.copy()
        if self.selection: colors[self.selection]=[1,.75,.2]
        data=np.column_stack((self.positions,radii,colors)).astype('f4')
        if self.atom_resources:
            buffer=self.atom_resources[1];buffer.orphan(data.nbytes);buffer.write(data.tobytes())
        else:
            buffer=self.ctx.buffer(data.tobytes())
            self.atom_vao=self.ctx.vertex_array(self.sphere,[(self.quad,'2f','corner'),(buffer,'3f 1f 3f /i','position','radius','color')])
            self.atom_resources=[self.atom_vao,buffer]
        if not self.bonds_dirty:
            self.dirty=False
            return
        pairs=self.molecule.bonds
        has_bonds=False
        if pairs is not None and len(pairs) and self.style!='Space filling':
            a,b=pairs[:,0].astype('i4'),pairs[:,1].astype('i4')
            # Zero-length bonds have no valid cylinder basis.
            valid=np.linalg.norm(self.positions[b]-self.positions[a],axis=1)>1e-5
            a,b=a[valid],b[valid]
            if len(a):
                has_bonds=True
                starts=self.positions[a];finishes=self.positions[b];ca=element_colors[a];cb=element_colors[b]
                orders=pairs[:,2][valid]
                # Additional instances share the same batch; camera motion changes only uniforms.
                extra=np.flatnonzero((orders==2)|(orders==3))
                if len(extra):
                    direction=finishes[extra]-starts[extra]
                    reference=np.tile(np.array([0,0,1],dtype='f4'),(len(extra),1))
                    parallel=np.abs(direction[:,2])/np.linalg.norm(direction,axis=1)>.9;reference[parallel]=[0,1,0]
                    offset=np.cross(direction,reference);offset/=np.linalg.norm(offset,axis=1,keepdims=True);offset*=.22
                    for multiplier in (-1,1):
                        starts=np.concatenate((starts,self.positions[a[extra]]+multiplier*offset))
                        finishes=np.concatenate((finishes,self.positions[b[extra]]+multiplier*offset))
                        ca=np.concatenate((ca,element_colors[a[extra]]));cb=np.concatenate((cb,element_colors[b[extra]]))
                    keep=np.ones(len(starts),dtype=bool);keep[extra[orders[extra]==2]]=False
                    starts,finishes,ca,cb=starts[keep],finishes[keep],ca[keep],cb[keep]
                aromatic=np.flatnonzero(orders==1.5)
                if len(aromatic):
                    direction=self.positions[b[aromatic]]-self.positions[a[aromatic]]
                    reference=np.tile(np.array([0,0,1],dtype='f4'),(len(aromatic),1))
                    reference[np.abs(direction[:,2])/np.linalg.norm(direction,axis=1)>.9]=[0,1,0]
                    offset=np.cross(direction,reference);offset/=np.linalg.norm(offset,axis=1,keepdims=True);offset*=.23
                    base=self.positions[a[aromatic]]+offset
                    for low,high in ((.08,.24),(.32,.48),(.56,.72),(.80,.96)):
                        starts=np.concatenate((starts,base+direction*low))
                        finishes=np.concatenate((finishes,base+direction*high))
                        color=element_colors[a[aromatic] if (low+high)/2<.5 else b[aromatic]]
                        ca=np.concatenate((ca,color));cb=np.concatenate((cb,color))
                data=np.column_stack((starts,finishes,ca,cb)).astype('f4')
                if self.bond_resources:
                    buffer=self.bond_resources[1];buffer.orphan(data.nbytes);buffer.write(data.tobytes())
                else:
                    buffer=self.ctx.buffer(data.tobytes())
                    self.bond_vao=self.ctx.vertex_array(self.bond,[(self.cylinder,'3f','vertex'),(buffer,'3f 3f 3f 3f /i','start','finish','color_a','color_b')])
                    self.coarse_bond_vao=self.ctx.vertex_array(self.bond,[(self.cylinder_coarse,'3f','vertex'),(buffer,'3f 3f 3f 3f /i','start','finish','color_a','color_b')])
                    self.bond_resources=[self.bond_vao,buffer,self.coarse_bond_vao]
                self.bond_count=len(starts)
        if not has_bonds:
            for r in self.bond_resources:r.release()
            self.bond_resources=[]
        self.dirty=False;self.bonds_dirty=False

    def _upload_surfaces(self):
        for r in self.surface_resources: r.release()
        self.surface_resources=[]; self.surface_vaos=[]
        for vertices,normals,color in self.surfaces:
            buffer=self.ctx.buffer(np.column_stack((vertices,normals)).astype('f4').tobytes())
            vao=self.ctx.vertex_array(self.surface,[(buffer,'3f 3f','position','normal_in')])
            self.surface_resources.extend([vao,buffer]); self.surface_vaos.append((vao,color))
        self.surface_dirty=False

    def paintGL(self):
        if self.ctx is None:
            painter=QPainter(self); painter.fillRect(self.rect(),QColor('#f8fafc'))
            painter.setPen(QColor('#9b2727')); painter.drawText(self.rect().adjusted(30,30,-30,-30),Qt.TextFlag.TextWordWrap,self.error or 'Initializing GPU…'); painter.end(); return
        t=time.perf_counter()
        fbo=self.ctx.detect_framebuffer(self.defaultFramebufferObject()); fbo.use()
        self.ctx.viewport=(0,0,round(self.width()*self.devicePixelRatioF()),round(self.height()*self.devicePixelRatioF()))
        self.ctx.clear(.965,.976,.988,1,depth=1)
        self.ctx.enable(moderngl.DEPTH_TEST); self.ctx.disable(moderngl.BLEND|moderngl.CULL_FACE)
        self.draw_calls=0
        if self.molecule is not None:
            if self.dirty: self._upload()
            if self.surface_dirty: self._upload_surfaces()
            rotation=self.rotation(); origin=self.origin.copy()
            origin+=rotation.T@np.array([*self.pan,0],dtype='f4')
            extent=self.extent()
            for program in (self.sphere,self.bond,self.surface):
                program['rotation'].write(rotation.T.tobytes()); program['origin'].value=tuple(origin); program['extent'].value=tuple(extent)
            bond_radius=.13 if self.style=='Sticks' else .105
            self.bond['bond_radius'].value=bond_radius
            if self.bond_resources:
                # Below a two-pixel radius, reduce cylinder tessellation. Keep
                # atom/bond identities, orders and colors; zooming restores detail.
                radius_pixels=bond_radius*self.height()*self.devicePixelRatioF()/(2*self.scale)
                coarse=radius_pixels<2 and not self.force_full_bonds
                vao=self.coarse_bond_vao if coarse else self.bond_vao
                vao.render(moderngl.TRIANGLES,instances=self.bond_count); self.draw_calls+=1
            self.atom_vao.render(moderngl.TRIANGLE_STRIP,vertices=4,instances=len(self.positions)); self.draw_calls+=1
            for vao,color in self.surface_vaos:
                self.surface['tint'].value=color; vao.render(moderngl.TRIANGLES); self.draw_calls+=1
        # Qt owns its default FBO. Do not call release() on the detected wrapper.
        self.paint_times.append((time.perf_counter()-t)*1000)
        self.paint_times=self.paint_times[-240:]

    def project(self):
        view=(self.positions-self.origin)@self.rotation().T; view[:,:2]-=self.pan
        ext=self.extent()
        xy=np.column_stack(((view[:,0]/ext[0]+1)*self.width()/2,(1-view[:,1]/ext[1])*self.height()/2))
        return xy,view[:,2]

    def mousePressEvent(self,event):
        self.setFocus(Qt.FocusReason.MouseFocusReason)
        self.last=event.position(); self.press=event.position(); self.dragged=False
        mode=self.effective_selection_mode()
        if event.button()==Qt.MouseButton.LeftButton and not event.modifiers() & Qt.KeyboardModifier.ShiftModifier and mode in ('rectangle','brush'):
            self.selection_drag=mode; self.brush_seen.clear()
            if mode=='rectangle':
                self.rubber_band.setGeometry(QRect(self.press.toPoint(),self.press.toPoint()))
                self.rubber_band.show()
            else:self.brush_at(event.position())

    def mouseMoveEvent(self,event):
        if self.last is None: return
        if self.selection_drag:
            if self.selection_drag=='rectangle':self.rubber_band.setGeometry(QRect(self.press.toPoint(),event.position().toPoint()).normalized())
            else:self.brush_at(event.position())
            return
        delta=event.position()-self.last; self.last=event.position()
        if (event.position()-self.press).manhattanLength()>4: self.dragged=True
        if not self.dragged:return
        if event.buttons() & Qt.MouseButton.MiddleButton or event.modifiers() & Qt.KeyboardModifier.ShiftModifier:
            self.pan+=np.array([-delta.x(),delta.y()],dtype='f4')*(2*self.scale/max(1,self.height()))
        elif event.buttons() & Qt.MouseButton.RightButton or event.modifiers() & Qt.KeyboardModifier.ControlModifier:
            self.roll+=delta.x()*.008
            self.scale=float(np.clip(self.scale*math.exp(delta.y()*.01),.25,100000))
        elif event.buttons() & Qt.MouseButton.LeftButton:
            self.yaw+=delta.x()*.008; self.pitch+=delta.y()*.008
        self.update()

    def mouseReleaseEvent(self,event):
        if self.selection_drag and event.button()==Qt.MouseButton.LeftButton:
            if self.selection_drag=='rectangle' and self.positions is not None:
                box=QRect(self.press.toPoint(),event.position().toPoint()).normalized()
                xy,_=self.project()
                inside=(xy[:,0]>=box.left())&(xy[:,0]<=box.right())&(xy[:,1]>=box.top())&(xy[:,1]<=box.bottom())
                indices=np.flatnonzero(inside).tolist()
                if box.width()<=4 and box.height()<=4:
                    hit=self.hit_atom(event.position()); indices=[hit] if hit>=0 else []
                self.selectionRequested.emit(indices,bool(event.modifiers() & Qt.KeyboardModifier.ControlModifier))
            self.cancel_selection_drag();return
        click=self.last is not None and not self.dragged and (event.position()-self.press).manhattanLength()<=4
        if click and self.positions is not None and event.button()==Qt.MouseButton.LeftButton:
            index=self.hit_atom(event.position())
            selecting=self.click_select_held or self.fragment_select_held
            if index>=0:
                if self.interaction_mode=='add' and not selecting:self.replacementRequested.emit(index)
                else:self.atomPicked.emit(index)
            elif self.interaction_mode=='add' and not selecting:
                _,z=self.project()
                depth=float(z[self.selection[0]]) if len(self.selection)==1 else 0.
                self.placementRequested.emit(self.screen_to_world(event.position(),depth))
            else:self.atomPicked.emit(-1)
        self.last=None

    def screen_to_world(self,point,depth=0.):
        """Place on a camera-facing plane, through the selected atom or view center."""
        ext=self.extent()
        view=np.array([(2*point.x()/max(1,self.width())-1)*ext[0]+self.pan[0],
                       (1-2*point.y()/max(1,self.height()))*ext[1]+self.pan[1],depth],dtype='f4')
        return (self.origin+self.rotation().T@view).astype('f4')

    def wheelEvent(self,event):
        self.scale=float(np.clip(self.scale*math.exp(-event.angleDelta().y()/120*.12),.25,100000)); self.update()

    def keyPressEvent(self,event):
        if event.key()==Qt.Key.Key_Home:self.fit()
        elif event.key()==Qt.Key.Key_F and not event.modifiers():self.fragment_select_held=True
        elif event.key()==Qt.Key.Key_Escape:
            self.cancel_selection_drag();self.held_selection_modes=[];self.fragment_select_held=False;self.click_select_held=False;self.exitPlacement.emit()
        elif event.key()==Qt.Key.Key_C and not event.modifiers():self.click_select_held=True
        elif event.key() in (Qt.Key.Key_R,Qt.Key.Key_B) and not event.modifiers() & (Qt.KeyboardModifier.ShiftModifier|Qt.KeyboardModifier.AltModifier|Qt.KeyboardModifier.MetaModifier):
            key=event.key();mode='rectangle' if key==Qt.Key.Key_R else 'brush'
            if not event.isAutoRepeat():self.held_selection_modes.append((key,mode))
        elif event.key()==Qt.Key.Key_A and not event.modifiers():self.selectAllRequested.emit()
        elif event.key()==Qt.Key.Key_N and not event.modifiers():self.clearSelectionRequested.emit()
        elif event.key() in (Qt.Key.Key_Delete,Qt.Key.Key_Backspace):self.deleteRequested.emit()
        elif event.key() in (Qt.Key.Key_Left,Qt.Key.Key_Right,Qt.Key.Key_Up,Qt.Key.Key_Down):
            dx={Qt.Key.Key_Left:-1,Qt.Key.Key_Right:1}.get(event.key(),0)
            dy={Qt.Key.Key_Up:-1,Qt.Key.Key_Down:1}.get(event.key(),0)
            if event.modifiers() & Qt.KeyboardModifier.ShiftModifier:
                self.pan+=np.array([-dx,dy],dtype='f4')*self.scale*.06
            elif event.modifiers() & Qt.KeyboardModifier.ControlModifier:
                self.roll+=dx*.08;self.scale=float(np.clip(self.scale*math.exp(dy*.08),.25,100000))
            else:self.yaw+=dx*.08;self.pitch+=dy*.08
            self.update()
        else: super().keyPressEvent(event)

    def keyReleaseEvent(self,event):
        if event.key()==Qt.Key.Key_C and not event.isAutoRepeat():self.click_select_held=False
        elif event.key()==Qt.Key.Key_F and not event.isAutoRepeat():self.fragment_select_held=False
        elif event.key() in (Qt.Key.Key_R,Qt.Key.Key_B) and not event.isAutoRepeat():
            self.held_selection_modes=[entry for entry in self.held_selection_modes if entry[0]!=event.key()]
        else:super().keyReleaseEvent(event)

    def focusOutEvent(self,event):
        self.click_select_held=False;self.fragment_select_held=False;self.held_selection_modes=[];self.cancel_selection_drag()
        super().focusOutEvent(event)
