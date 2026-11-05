"""Element coordination templates and transactional fragment placement.

The catalog is derived from the supplied GaussView elements.eftb and .frg
resources. Six carbon patterns use independently constructed ideal directions;
other patterns use normalized resource directions and equilibrium bond lengths.
Initial placement does not reproduce the original molecular mechanics Clean.
Zero denotes an editable open valence, never element 114 or a real atom.
"""
from dataclasses import dataclass
from itertools import permutations
from importlib.resources import files
import json
import math
import numpy as np
from .model import RADII, SYMBOLS


@dataclass(frozen=True)
class CarbonFragment:
    key: str
    name: str
    orders: tuple
    directions: tuple


T = np.sqrt(3)/2
CARBON_FRAGMENTS = (
    CarbonFragment('atom', 'Carbon atom', (), ()),
    CarbonFragment('st', 'Carbon · linear (single–triple)', (1, 3), ((-1,0,0),(1,0,0))),
    CarbonFragment('dd', 'Carbon · linear (double–double)', (2, 2), ((-1,0,0),(1,0,0))),
    CarbonFragment('ssd', 'Carbon · trigonal (single–single–double)', (1, 1, 2), ((-1,0,0),(.5,-T,0),(.5,T,0))),
    CarbonFragment('saa', 'Carbon · aromatic (single–aromatic–aromatic)', (1, 1.5, 1.5), ((-1,0,0),(.5,-T,0),(.5,T,0))),
    CarbonFragment('tet', 'Carbon · tetrahedral', (1,1,1,1), ((-1,0,0),(1/3,np.sqrt(8/9),0),(1/3,-np.sqrt(2/9),np.sqrt(2/3)),(1/3,-np.sqrt(2/9),-np.sqrt(2/3)))),
)
FRAGMENTS = {f.key:f for f in CARBON_FRAGMENTS}
CATALOG = json.loads(files('openview').joinpath('data/element_fragments.json').read_text())
for _template in CATALOG['templates']:
    if _template['key'] not in FRAGMENTS:
        FRAGMENTS[_template['key']] = CarbonFragment(_template['key'],_template['name'],tuple(_template['orders']),tuple(map(tuple,_template['directions'])))


def element_options(number):
    return CATALOG['elements'].get(SYMBOLS[number],{'name':SYMBOLS[number],'default':'atom','templates':['atom']})


def fragment_name(number,key):
    template=next(t for t in CATALOG['templates'] if t['key']==key)
    return element_options(number)['name']+' · '+template['name']


def available_fragments(number):
    return tuple(FRAGMENTS[key] for key in element_options(number)['templates'])


def unit(v):
    length=np.linalg.norm(v)
    if length<1e-8: raise ValueError('Cannot attach at coincident atoms.')
    return np.asarray(v,dtype=float)/length


def align(source,target):
    """Proper rotation between two directions, including antiparallel ones."""
    a,b=unit(source),unit(target); v=np.cross(a,b); c=float(a@b)
    if c < -1+1e-8:
        axis=unit(np.cross(a,np.eye(3)[np.argmin(abs(a))]))
        return 2*np.outer(axis,axis)-np.eye(3)
    k=np.array([[0,-v[2],v[1]],[v[2],0,-v[0]],[-v[1],v[0],0]])
    return np.eye(3)+k+k@k/(1+c)


def bond_length(other,order,center=6):
    if 0 in (other,center):return 1.0
    radii=CATALOG['covalent_radii']
    return float(radii.get(SYMBOLS[other],RADII[other])+radii.get(SYMBOLS[center],RADII[center])) * {1:1.,1.5:.91,2:.88,3:.78}[order]


def axial_rotation(axis,angle):
    x,y,z=unit(axis);k=np.array([[0,-z,y],[z,0,-x],[-y,x,0]])
    return np.eye(3)+np.sin(angle)*k+(1-np.cos(angle))*(k@k)


def orient_join(rotation,directions,used,parent,order,center,positions,bonds,drop):
    """Resolve the free torsion using molecular geometry, not camera roll."""
    free=[i for i in range(len(directions)) if i!=used]
    if not free:return rotation
    axis=unit(positions[parent]-center)
    local=directions[free]@rotation.T
    across=local-np.outer(local@axis,axis)
    if np.max(np.linalg.norm(across,axis=1))<1e-8:return rotation
    if order in (1.5,2):
        neighbors=[int(b if int(a)==parent else a) for a,b,_ in bonds if parent in (int(a),int(b)) and int(b if int(a)==parent else a) not in drop]
        if neighbors:
            reference=positions[neighbors]-positions[parent]
            reference-=np.outer(reference@axis,axis)
            lengths=np.linalg.norm(reference,axis=1)
            if lengths.max()>1e-8:
                a=unit(across[np.argmax(np.linalg.norm(across,axis=1))]);b=unit(reference[np.argmax(lengths)])
                angle=np.arctan2(axis@np.cross(a,b),a@b)
                return axial_rotation(axis,angle)@rotation
    # For single bonds, choose the least crowded of 24 torsions. This is a
    # geometric starting placement, not the original molecular mechanics Clean.
    nearby=[i for i in range(len(positions)) if i not in drop and i!=parent and np.linalg.norm(positions[i]-center)<5]
    if not nearby:return rotation
    best_score=-np.inf;best=rotation
    for angle in np.linspace(0,2*np.pi,24,endpoint=False):
        r=axial_rotation(axis,angle)@rotation
        trial=center+(directions[free]@r.T)*1.09
        distances=np.sum((trial[:,None,:]-positions[nearby][None,:,:])**2,axis=2)
        score=float(np.min(distances))
        if score>best_score+1e-8:best_score=score;best=r
    return best


def build_element(numbers,frames,bonds,key,*,atomic_number,position=None,target=None,orientation=None,match_orders=False):
    """Return new topology and central index; input arrays are never mutated.

    Empty-space placement creates a distinct fragment. A clicked terminal H
    or open valence is consumed to attach the fragment. A clicked heavy atom
    keeps its nonterminal/heavy neighbors, replacing attached H/open valences.
    Target bond orders are retained (original built-in default). Strict join
    matching can be requested separately. Ring and
    existing heavy-atom coordinates remain intact during center replacement.
    """
    fragment=FRAGMENTS[key]
    if key not in element_options(atomic_number)['templates']:raise ValueError('This coordination pattern is not available for the selected element.')
    numbers=np.asarray(numbers,dtype='i4'); frames=np.asarray(frames,dtype='f4')
    bonds=np.asarray(bonds,dtype='f4').reshape(-1,3)
    n=len(numbers); rotation=np.eye(3) if orientation is None else np.asarray(orientation,dtype=float)
    if rotation.shape!=(3,3) or not np.isfinite(rotation).all(): raise ValueError('Invalid placement orientation.')
    degrees=np.bincount(bonds[:,:2].astype(int).ravel(),minlength=n)
    drop=set(); attached=[]; replacing=False
    if target is None:
        center=np.asarray(position,dtype=float)
        if center.shape!=(3,) or not np.isfinite(center).all(): raise ValueError('Invalid placement position.')
        centers=np.tile(center,(len(frames),1))
    else:
        if not 0<=target<n: raise ValueError('Unknown target atom.')
        edges=[(int(b if int(a)==target else a),float(order)) for a,b,order in bonds if target in (int(a),int(b))]
        if numbers[target] in (0,1) and len(edges)==1:
            parent,order=edges[0]
            if numbers[parent]==0: raise ValueError('Attach to the open valence of a real atom.')
            direction=np.array([unit(p[target]-p[parent]) for p in frames])
            centers=frames[:,parent,:]+direction*bond_length(int(numbers[parent]),order,atomic_number)
            drop.add(target); attached=[(parent,order)]
        elif numbers[target]==0:
            raise ValueError('An open valence must have exactly one bond.')
        else:
            replacing=True; centers=frames[:,target,:].copy()
            for neighbor,order in edges:
                if numbers[neighbor] in (0,1) and degrees[neighbor]==1: drop.add(neighbor)
                else: attached.append((neighbor,order))
    orders=fragment.orders
    # Bare-atom replacement intentionally keeps all old bonds/terminal atoms.
    if not orders and replacing:
        result=numbers.copy(); result[target]=atomic_number
        return result,frames.copy(),bonds.copy(),target
    # Enumerating ten-coordinate metal permutations would freeze the UI.
    # Use bounded assignment search for large coordination numbers instead.
    if len(attached)>len(orders):choices=[]
    elif math.perm(len(orders),len(attached))<=1000:
        choices=list(permutations(range(len(orders)),len(attached)))
    else:
        desired=np.array([unit(frames[0,a]-centers[0]) for a,_ in attached])
        directions=np.asarray(fragment.directions)
        candidates=[]
        for first in range(len(orders)):
            r=align(directions[first],desired[0]);oriented=directions@r.T
            costs=(100*np.abs(np.asarray([order for _,order in attached])[:,None]-np.asarray(orders)[None,:])
                   +np.sum((desired[:,None,:]-oriented[None,:,:])**2,axis=2)).tolist()
            # Dynamic programming over occupied ports: O(k*n*2**n).
            states={0:(0.,())}
            for j,(_,order) in enumerate(attached):
                next_states={}
                for mask,(cost,path) in states.items():
                    for port in range(len(orders)):
                        if mask&(1<<port):continue
                        score=cost+costs[j][port]
                        new_mask=mask|(1<<port)
                        if new_mask not in next_states or score<next_states[new_mask][0]:next_states[new_mask]=(score,path+(port,))
                states=next_states
            candidates.append(min(states.values())[1])
        choices=list(dict.fromkeys(candidates))
    if choices and attached:
        mismatches=[sum(abs(orders[port]-order) for port,(_,order) in zip(p,attached)) for p in choices]
        best_match=min(mismatches)
        choices=[p for p,score in zip(choices,mismatches) if score==best_match and (not match_orders or score<1e-6)]
    if attached and orders and not choices:
        raise ValueError('This element pattern has no matching free bond. Choose a pattern with sufficient coordination and matching bond types if strict matching is enabled.')
    if not attached: choices=[()]
    if not orders: choices=[()]  # a bare atom may terminate any existing bond
    directions=np.asarray(fragment.directions,dtype=float).reshape(-1,3)
    mapping=choices[0]
    rotations=[]
    for frame_index,p in enumerate(frames):
        if attached and orders:
            desired=np.array([unit(p[a]-centers[frame_index]) for a,_ in attached])
            best=None
            for candidate in choices if frame_index==0 else [mapping]:
                source=directions[list(candidate)]
                if len(candidate)==1:
                    fitted=align(rotation@source[0],desired[0])@rotation
                else:
                    u,_,vt=np.linalg.svd(source.T@desired)
                    correction=np.eye(3);correction[-1,-1]=np.linalg.det(u@vt)
                    fitted=(u@correction@vt).T
                error=float(np.sum((source@fitted.T-desired)**2))
                if best is None or error<best[0]: best=(error,candidate,fitted)
            _,mapping,fitted=best
            if len(attached)==1 and not replacing:
                parent,order=attached[0]
                fitted=orient_join(fitted,directions,mapping[0],parent,order,centers[frame_index],p,bonds,drop)
            rotations.append(fitted)
        else: rotations.append(rotation)
    keep=np.array([i for i in range(n) if i not in drop],dtype=int)
    remap=np.full(n,-1,dtype=int); remap[keep]=np.arange(len(keep))
    result_numbers=list(numbers[keep]); result_frames=frames[:,keep,:].copy()
    if replacing: center_index=int(remap[target]);result_numbers[center_index]=atomic_number
    else:
        center_index=len(result_numbers);result_numbers.append(atomic_number)
        result_frames=np.concatenate((result_frames,centers[:,None,:]),axis=1)
    result_bonds=[]
    for a,b,order in bonds:
        a,b=int(a),int(b)
        if a not in drop and b not in drop: result_bonds.append([remap[a],remap[b],order])
    if not replacing:
        for neighbor,order in attached: result_bonds.append([remap[neighbor],center_index,order])
    for port,order in enumerate(orders):
        if port in mapping: continue
        index=len(result_numbers);result_numbers.append(1 if order==1 else 0)
        # Open-valence endpoint is a visible construction handle, not an atom
        # with a chemical bond length. Finished bonds get their own length.
        length=bond_length(1,1,atomic_number) if order==1 else 1.0
        xyz=np.array([c+r@directions[port]*length for c,r in zip(centers,rotations)])
        result_frames=np.concatenate((result_frames,xyz[:,None,:]),axis=1)
        result_bonds.append([center_index,index,order])
    return (np.asarray(result_numbers,dtype='i4'),result_frames.astype('f4'),
            np.asarray(result_bonds,dtype='f4').reshape(-1,3),center_index)


def build_carbon(numbers,frames,bonds,key,**kwargs):
    """Compatibility wrapper for existing carbon builder callers."""
    return build_element(numbers,frames,bonds,key,atomic_number=6,**kwargs)
