from dataclasses import dataclass, field
from itertools import product
import numpy as np

SYMBOLS = "X H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og".split()
ATOMIC_NUMBER = {s.lower(): i for i, s in enumerate(SYMBOLS)}
RADII = np.full(119, 1.4, dtype='f4')
RADII[0] = .3  # internal open-valence construction handle
for z,r in {1:.31,2:.28,3:1.28,4:.96,5:.84,6:.76,7:.71,8:.66,9:.57,10:.58,11:1.66,12:1.41,13:1.21,14:1.11,15:1.07,16:1.05,17:1.02,18:1.06,19:2.03,20:1.76,26:1.32,29:1.32,30:1.22,35:1.20,53:1.39}.items(): RADII[z]=r
COLORS = np.tile(np.array([.55,.62,.72],dtype='f4'), (119,1))
COLORS[0] = [1., .48, .12]
for z,c in {1:(.92,.94,1),6:(.35,.4,.48),7:(.23,.45,1),8:(1,.24,.3),9:(.3,.85,.5),15:(1,.6,.2),16:(1,.84,.2),17:(.23,.85,.4),26:(.85,.5,.25),35:(.66,.26,.2),53:(.6,.3,.83)}.items(): COLORS[z]=c

def element(value):
    s=str(value).strip()
    if s.isdigit() and 1 <= int(s) < len(SYMBOLS): return int(s)
    z=ATOMIC_NUMBER.get(s.lower())
    if not z: raise ValueError(f"Unknown element: {value!r}")
    return z

def connected_fragment(bonds, atom_count, start):
    """Return a connected component in atom order, including isolated atoms."""
    if not 0 <= start < atom_count:
        return []
    adjacency = [[] for _ in range(atom_count)]
    for a, b, order in bonds:
        if order > 0:
            a, b = int(a), int(b)
            adjacency[a].append(b); adjacency[b].append(a)
    seen = {start}; pending = [start]
    while pending:
        for other in adjacency[pending.pop()]:
            if other not in seen:
                seen.add(other); pending.append(other)
    return sorted(seen)

@dataclass
class Molecule:
    numbers: np.ndarray
    frames: np.ndarray
    title: str = "Untitled"
    bonds: np.ndarray | None = None
    energies: list = field(default_factory=list)
    frequencies: list = field(default_factory=list)
    intensities: list = field(default_factory=list)
    modes: np.ndarray | None = None
    charge: int = 0
    multiplicity: int = 1
    metadata: dict = field(default_factory=dict)
    grid: object = None

    def __post_init__(self):
        self.numbers=np.asarray(self.numbers,dtype='i4')
        self.frames=np.asarray(self.frames,dtype='f4')
        if self.frames.ndim==2: self.frames=self.frames[None,:,:]
        if not len(self.numbers) or self.frames.ndim!=3 or self.frames.shape[1:]!=(len(self.numbers),3):
            raise ValueError("Coordinates must contain one xyz position per atom in every frame")
        if not np.isfinite(self.frames).all(): raise ValueError("Non-finite coordinates")
        if np.any(self.numbers<0) or np.any(self.numbers>=len(SYMBOLS)): raise ValueError("Invalid atomic number")
        if self.bonds is not None:
            self.bonds=np.asarray(self.bonds,dtype='f4').reshape(-1,3)
            if not np.isfinite(self.bonds).all() or np.any(self.bonds[:,:2]!=np.floor(self.bonds[:,:2])): raise ValueError('Invalid bond indices/orders')
            if len(self.bonds) and (self.bonds[:,:2].min()<0 or self.bonds[:,:2].max()>=len(self.numbers)): raise ValueError("Invalid bond index")

    def infer_bonds(self):
        if self.bonds is None: self.bonds=infer_bonds(self.numbers,self.frames[0])
        return self

def infer_bonds(numbers, positions, tolerance=1.2):
    """Spatial hash: bounded-neighborhood search, no NxN distance matrix."""
    n=len(numbers)
    if not n: return np.empty((0,3),dtype='i4')
    radii=RADII[numbers]; cell=float(2*radii.max()*tolerance)
    keys=np.floor(positions/cell).astype('i8')
    buckets={}; bonds=[]
    offsets=list(product((-1,0,1),repeat=3))
    for i,key in enumerate(keys):
        candidates=[]
        for dx,dy,dz in offsets: candidates.extend(buckets.get((key[0]+dx,key[1]+dy,key[2]+dz),()))
        if candidates:
            js=np.asarray(candidates)
            d2=np.sum((positions[js]-positions[i])**2,axis=1)
            hit=js[(d2>.16)&(d2<(tolerance*(radii[js]+radii[i]))**2)]
            bonds.extend((int(j),i,1) for j in hit)
        buckets.setdefault(tuple(key),[]).append(i)
    return np.asarray(bonds,dtype='i4').reshape(-1,3)

def measure(positions, indices):
    if len(indices)>4:return f'{len(indices):,} atoms selected · select 2–4 to measure'
    p=np.asarray(positions,dtype='f8')[indices]
    if len(p)==2: return f"Distance  {np.linalg.norm(p[1]-p[0]):.4f} Å"
    if len(p)==3:
        a,b=p[0]-p[1],p[2]-p[1]; denom=np.linalg.norm(a)*np.linalg.norm(b)
        if denom<1e-12: return "Angle undefined (coincident atoms)"
        return f"Angle  {np.degrees(np.arccos(np.clip(a@b/denom,-1,1))):.3f}°"
    if len(p)==4:
        b0=-(p[1]-p[0]); b1=p[2]-p[1]; b2=p[3]-p[2]
        if np.linalg.norm(b1)<1e-12: return "Dihedral undefined"
        b1/=np.linalg.norm(b1); v=b0-(b0@b1)*b1; w=b2-(b2@b1)*b1
        if min(np.linalg.norm(v),np.linalg.norm(w))<1e-12: return "Dihedral undefined (collinear atoms)"
        return f"Dihedral  {np.degrees(np.arctan2(np.cross(b1,v)@w,v@w)):.3f}°"
    return "Select 2–4 atoms to measure"
