"""Local geometry q and reference-projected lambda observables."""
from __future__ import annotations
import numpy as np
from ..io.trajectory import Frame
from ..structure.neighbors import dynamic_cation_n_neighbors

def compute_q(frame: Frame,cutoff_A:float=3.5):
    """Dynamic four-N tetrahedral c-axis third moment."""
    cation,_n,vectors=dynamic_cation_n_neighbors(frame,4,cutoff_A)
    axis=frame.lattice[2]/np.linalg.norm(frame.lattice[2]); proj=np.einsum('ijk,k->ij',vectors,axis); cube=proj**3; den=np.abs(cube).sum(axis=1)
    q=np.divide(cube.sum(axis=1),den,out=np.zeros_like(den),where=den>1e-14)
    return cation[:,0],q

def compute_lambda(frame:Frame,reference:dict):
    """Project MIC fractional displacements onto a matched endpoint path.

    Accepts either `frac0`/`disp` from the existing N0/N1 reference.npz, or
    `parent_positions`/`reversed_positions`/`cell` Cartesian endpoint arrays.
    Atom ordering must match the trajectory.
    """
    current=frame.positions@np.linalg.inv(frame.lattice)
    if 'frac0' in reference and 'disp' in reference:
        frac0=np.asarray(reference['frac0'],float);path=np.asarray(reference['disp'],float)
        if frac0.shape!=current.shape or path.shape!=current.shape:raise ValueError('fractional reference atom count/order must match trajectory')
    else:
        parent=np.asarray(reference['parent_positions'],float);rev=np.asarray(reference['reversed_positions'],float);cell=np.asarray(reference['cell'],float).reshape(3,3)
        if parent.shape!=frame.positions.shape or rev.shape!=frame.positions.shape:raise ValueError('reference endpoint atom counts/order must match trajectory')
        frac0=parent@np.linalg.inv(cell);target=rev@np.linalg.inv(cell);path=target-frac0;path-=np.rint(path)
    delta=current-frac0;delta-=np.rint(delta)
    den=np.einsum('ij,ij->i',path,path)
    return np.divide(np.einsum('ij,ij->i',delta,path),den,out=np.zeros_like(den),where=den>1e-12)

def classify_q(q:np.ndarray,threshold:float=.2):
    phase=np.zeros(np.asarray(q).shape,np.int8); phase[np.asarray(q)>threshold]=1; phase[np.asarray(q)<-threshold]=-1; return phase
def classify_lambda(value:np.ndarray,threshold:float=.5):
    """The target endpoint side is +1; parent side is -1."""
    return np.where(np.asarray(value)>threshold,1,-1).astype(np.int8)
def classify_lambda_grid(frame:Frame,value:np.ndarray,threshold:float=.5,shape=(24,24)):
    """Average atomic lambda on periodic in-plane cells then classify cells.

    This reproduces the established N0/N1 observable: every atom contributes
    to its fractional a-b cell, and a cell is reversed when mean lambda > 0.5.
    Returned phase labels are a 2-D grid, not atom counts.
    """
    nx,ny=map(int,shape);frac=frame.positions@np.linalg.inv(frame.lattice);frac-=np.floor(frac)
    ix=np.minimum((frac[:,0]*nx).astype(int),nx-1);iy=np.minimum((frac[:,1]*ny).astype(int),ny-1)
    flat=ix*ny+iy;total=nx*ny;sums=np.bincount(flat,weights=np.asarray(value),minlength=total);counts=np.bincount(flat,minlength=total)
    mean=np.divide(sums,counts,out=np.zeros(total),where=counts>0);phase=np.where(mean>threshold,1,-1).astype(np.int8)
    return phase.reshape(nx,ny),mean.reshape(nx,ny),counts.reshape(nx,ny)

def phase_fractions(phase:np.ndarray):
    phase=np.asarray(phase)
    if not phase.size:return dict(plus_fraction=float('nan'),minus_fraction=float('nan'),transition_fraction=float('nan'))
    return dict(plus_fraction=float(np.mean(phase==1)),minus_fraction=float(np.mean(phase==-1)),transition_fraction=float(np.mean(phase==0)))
