"""GPUMD thermo.out parsing and rolling robust anomaly detection."""
from __future__ import annotations
import numpy as np
DEFAULT_COLUMNS={'temperature':0,'kinetic_energy':1,'potential_energy':2,'sigma_zz':8,'c_A':18}
def read_thermo(path,columns=None):
    cols={**DEFAULT_COLUMNS,**(columns or {})};rows=[]
    with open(path,'rt',encoding='utf-8') as f:
        for line in f:
            vals=np.fromstring(line,sep=' ')
            if vals.size:rows.append(vals)
    if not rows:return {'columns':cols,'values':np.empty((0,0))}
    a=np.full((len(rows),max(map(len,rows))),np.nan)
    for i,row in enumerate(rows):a[i,:len(row)]=row
    return {'columns':cols,'values':a}
def check_thermo(table,interval_fs=2.5,max_temperature_K=2000.,max_pe_jump_eV_atom=.5,natoms=32768):
    a=table['values'];cols=table['columns'];reasons=[];warnings=[];first=None
    if not len(a):return {'status':'WARNING','reasons':['thermo.out contains no numeric rows'],'warnings':[],'first_bad_time_ps':None,'metrics':{}}
    def col(name):
        ix=cols.get(name);return a[:,ix] if ix is not None and ix<a.shape[1] else np.full(len(a),np.nan)
    temp,pe,sig=col('temperature'),col('potential_energy'),col('sigma_zz')
    def fail(i,reason):
        nonlocal first
        if first is None:first={'row':int(i),'time_ps':float(i*interval_fs/1000.),'reason':reason}
        if reason not in reasons:reasons.append(reason)
    bad=np.flatnonzero(~np.isfinite(temp))
    if bad.size:fail(int(bad[0]),'non-finite temperature')
    for i in np.flatnonzero(np.isfinite(temp)):
        hist=temp[max(0,int(i)-101):i];hist=hist[np.isfinite(hist)]
        if temp[i]>max_temperature_K:fail(int(i),'temperature exceeds configured absolute limit')
        if hist.size>=10:
            med=float(np.median(hist));mad=float(np.median(np.abs(hist-med)));scale=max(1.4826*mad,15.)
            if temp[i]>med+6*scale and temp[i]>800.:fail(int(i),'temperature runaway (rolling median/MAD)')
    jumps=np.flatnonzero(np.isfinite(pe[1:])&np.isfinite(pe[:-1])&(np.abs(np.diff(pe))/max(natoms,1)>max_pe_jump_eV_atom))
    if jumps.size:fail(int(jumps[0]+1),'potential-energy discontinuity')
    for i in np.flatnonzero(np.isfinite(sig)):
        hist=sig[max(0,int(i)-101):i];hist=hist[np.isfinite(hist)]
        if hist.size>=10:
            med=float(np.median(hist));mad=float(np.median(np.abs(hist-med)));scale=max(1.4826*mad,1e-8)
            if abs(sig[i]-med)>8*scale:
                warnings.append('sigma_zz rolling median/MAD outlier');break
    if len(a)<10:warnings.append('fewer than 10 thermo samples')
    return {'status':'FAILED' if reasons else 'WARNING' if warnings else 'PASS','reasons':reasons,'warnings':warnings,'first_bad_time_ps':first['time_ps'] if first else None,'first_bad':first,
            'metrics':{'n_samples':int(len(a)),'temperature_min_K':float(np.nanmin(temp)),'temperature_max_K':float(np.nanmax(temp)),'potential_energy_initial_eV':float(pe[0]),'potential_energy_final_eV':float(pe[-1]),'sigma_zz_min':float(np.nanmin(sig)),'sigma_zz_max':float(np.nanmax(sig))}}
