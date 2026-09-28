"""Dynamic nearest-N coordination in triclinic periodic cells via ASE."""
from __future__ import annotations
import numpy as np
from ase import Atoms
from ase.neighborlist import neighbor_list
from ..io.trajectory import Frame
CODE_TO_SYMBOL=np.asarray(['X','Al','Sc','N'])
def _atoms(frame: Frame) -> Atoms:
    return Atoms(CODE_TO_SYMBOL[frame.species],positions=frame.positions,cell=frame.lattice,pbc=True)
def dynamic_cation_n_neighbors(frame: Frame, coordination: int=4, cutoff_A: float=3.5):
    """Return cation and four nearest-N indices plus minimum-image vectors."""
    atoms=_atoms(frame); i,j,d=neighbor_list('ijD',atoms,cutoff_A,self_interaction=False)
    mask=((frame.species[i]==1)|(frame.species[i]==2))&(frame.species[j]==3)
    i,j,d=i[mask],j[mask],d[mask]
    order=np.lexsort((np.linalg.norm(d,axis=1),i)); i,j,d=i[order],j[order],d[order]
    cations=np.flatnonzero((frame.species==1)|(frame.species==2)); counts=np.bincount(i,minlength=frame.natoms)
    if np.any(counts[cations]<coordination):
        atom=int(cations[np.flatnonzero(counts[cations]<coordination)[0]])
        raise ValueError(f'fewer than {coordination} N neighbours within {cutoff_A:g} A for atom {atom}')
    take=np.zeros(len(i),bool); seen=np.zeros(frame.natoms,np.int8)
    for k,atom in enumerate(i):
        if seen[atom]<coordination: take[k]=True; seen[atom]+=1
    return i[take].reshape(-1,coordination),j[take].reshape(-1,coordination),d[take].reshape(-1,coordination,3)
def minimum_distances(frame: Frame, cutoff_A: float=4.0) -> dict[str,float]:
    atoms=_atoms(frame); i,j,d=neighbor_list('ijD',atoms,cutoff_A,self_interaction=False); dist=np.linalg.norm(d,axis=1)
    si,sj=frame.species[i],frame.species[j]
    masks={'min_AlN_A':((si==1)&(sj==3))|((si==3)&(sj==1)),
           'min_ScN_A':((si==2)&(sj==3))|((si==3)&(sj==2)),
           'min_cation_cation_A':((si==1)|(si==2))&((sj==1)|(sj==2))&(i<j)}
    return {key:(float(dist[mask].min()) if np.any(mask) else float('nan')) for key,mask in masks.items()}
