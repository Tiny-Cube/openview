"""Cartesian construction from Gaussian internal coordinates (no expression eval)."""
import math
import re
import numpy as np
from .model import element

def build_zmatrix(rows,tail):
    variables={};labels={};positions=[];numbers=[]
    for line in tail:
        match=re.fullmatch(r'\s*([A-Za-z]\w*)\s*(?:=|\s)\s*([-+\d.eEdD]+)\s*',line)
        if match:variables[match[1].lower()]=float(match[2].replace('D','E').replace('d','e'))
    def value(token):
        try:return float(token.replace('D','E').replace('d','e'))
        except ValueError:
            sign=-1 if token.startswith('-') else 1;key=token.lstrip('+-').lower()
            if key not in variables:raise ValueError(f'Undefined Z-matrix variable: {token}')
            return sign*variables[key]
    def reference(token,index):
        if token.isdigit():result=int(token)-1
        else:result=labels.get(token.lower(),-1)
        if result is None or not 0<=result<index:raise ValueError(f'Invalid Z-matrix reference: {token}')
        return result
    for index,row in enumerate(rows):
        expected=1+2*min(index,3)
        if len(row)!=expected:raise ValueError(f'Unsupported Z-matrix row {index+1}; expected {expected} fields')
        label=row[0];match=re.fullmatch(r'([A-Za-z]{1,2})\d*',label)
        if not match:raise ValueError(f'Invalid Z-matrix element: {label}')
        symbol=match[1];z=0 if symbol.lower()=='x' else element(symbol)
        if not index:pos=np.zeros(3)
        else:
            a=reference(row[1],index);r=value(row[2])
            if not np.isfinite(r) or r<=0:raise ValueError('Z-matrix distances must be positive')
            if index==1:pos=positions[a]+np.array([0,0,r])
            else:
                b=reference(row[3],index);theta=math.radians(value(row[4]))
                if a==b or not 0<=theta<=math.pi:raise ValueError('Invalid Z-matrix angle/reference')
                axis=positions[b]-positions[a];length=np.linalg.norm(axis)
                if length<1e-12:raise ValueError('Coincident Z-matrix references')
                axis/=length
                if index==2:
                    aux=np.array([1.,0,0]) if abs(axis[0])<.9 else np.array([0.,1,0]);normal=np.cross(axis,aux);phi=0
                else:
                    c=reference(row[5],index)
                    if len({a,b,c})<3:raise ValueError('Repeated Z-matrix references')
                    normal=np.cross(axis,positions[c]-positions[b]);phi=math.radians(value(row[6]))
                length=np.linalg.norm(normal)
                if length<1e-12:raise ValueError('Collinear Z-matrix references; add a dummy center')
                normal/=length;side=np.cross(normal,axis)
                pos=positions[a]+r*(math.cos(theta)*axis+math.sin(theta)*(math.cos(phi)*side-math.sin(phi)*normal))
        labels[label.lower()]=None if label.lower() in labels else index
        positions.append(pos);numbers.append(z)
    mask=np.asarray(numbers)>0
    return np.asarray(numbers)[mask],np.asarray(positions)[mask]
