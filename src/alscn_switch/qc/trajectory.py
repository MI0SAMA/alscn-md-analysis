"""Sampled geometric checks for trajectory frames."""
import numpy as np
from ..structure.neighbors import minimum_distances
def check_frame(frame,previous_c=None,max_c_jump_A=.5,check_distances=True):
    reasons=[];warnings=[];c=float(np.linalg.norm(frame.lattice[2]))
    if not np.all(np.isfinite(frame.positions)):reasons.append('non-finite position')
    if previous_c is not None and abs(c-previous_c)>max_c_jump_A:reasons.append('c-axis cell jump')
    good=np.isfinite(frame.forces).all(axis=1);fmax=float(np.linalg.norm(frame.forces[good],axis=1).max()) if good.any() else float('nan')
    if not good.all():warnings.append('force property unavailable or non-finite for some atoms')
    goodbec=np.isfinite(frame.bec);becmax=float(np.abs(frame.bec[goodbec]).max()) if goodbec.any() else float('nan')
    d=minimum_distances(frame) if check_distances else {}
    if np.isfinite(d.get('min_AlN_A',np.nan)) and d['min_AlN_A']<1.2:reasons.append('Al-N short contact')
    if np.isfinite(d.get('min_ScN_A',np.nan)) and d['min_ScN_A']<1.2:reasons.append('Sc-N short contact')
    return {'frame':frame.index,'time_ps':frame.time_fs/1000.,'c_A':c,'fmax_eV_A':fmax,'bec_abs_max':becmax,**d,'reasons':reasons,'warnings':warnings}
