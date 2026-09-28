"""Single-run streaming analysis and provenance outputs."""
from __future__ import annotations
import csv,json,os,platform,time
from pathlib import Path
import numpy as np
from . import __version__
from .io.manifest import load_run,load_yaml
from .io.trajectory import iter_extxyz
from .polarization.order import compute_q,compute_lambda,classify_q,classify_lambda_grid,phase_fractions
from .qc.thermo import read_thermo,check_thermo
from .qc.trajectory import check_frame

def _load_reference(path):
    if not path:raise ValueError('lambda observable requires analysis.reference NPZ file')
    with np.load(path,allow_pickle=False) as z:return {k:z[k] for k in z.files}

def analyze_run(manifest_path,outdir,backend='auto',stride=None,allow_failed=False):
    cfg=load_run(manifest_path); d=cfg.data; out=Path(outdir);out.mkdir(parents=True,exist_ok=True)
    analysis=d['analysis'];sampling=int(stride or analysis.get('frame_stride',1));observable=analysis.get('observable','q')
    if observable not in ('q','lambda'):raise ValueError('analysis.observable must be q or lambda')
    threshold=float(analysis.get('threshold',.2 if observable=='q' else .5));reference=_load_reference(analysis.get('reference')) if observable=='lambda' else None
    table=None
    thermo_qc={'status':'WARNING','reasons':['thermo file not provided'],'warnings':[],'first_bad_time_ps':None,'metrics':{}}
    thermo_path=d['files'].get('thermo')
    if thermo_path and Path(thermo_path).exists():
        table=read_thermo(thermo_path,d.get('qc',{}).get('thermo_columns'))
        thermo_qc=check_thermo(table,float(d.get('qc',{}).get('thermo_interval_fs',2.5)),float(d.get('qc',{}).get('max_temperature_K',2000)),float(d.get('qc',{}).get('max_pe_jump_eV_atom',.5)),int(d.get('structure',{}).get('natoms',32768)))
    fields=['run_id','frame','time_ps','observable','threshold','plus_fraction','minus_fraction','transition_fraction','fraction_basis','mean_order','max_order']
    tmp=out/'timeseries.csv.partial';qc_rows=[];previous_c=None;started=time.perf_counter();nread=0
    with tmp.open('w',newline='',encoding='utf-8') as f:
        writer=csv.DictWriter(f,fieldnames=fields);writer.writeheader()
        for frame in iter_extxyz(d['files']['trajectory'],sampling,backend):
            nread+=1;item=check_frame(frame,previous_c,float(analysis.get('max_c_jump_A',.5)),bool(analysis.get('check_distances',True)));previous_c=item['c_A'];qc_rows.append(item)
            if thermo_qc['status']=='FAILED':continue
            if observable=='q':_idx,val=compute_q(frame,float(analysis.get('neighbor_cutoff_A',3.5)));phase=classify_q(val,threshold)
            else:
                val=compute_lambda(frame,reference);phase,grid_lambda,_occupancy=classify_lambda_grid(frame,val,threshold,analysis.get('phase_grid',(24,24)))
            fractions=phase_fractions(phase)
            basis='cation' if observable=='q' else 'periodic_grid_cell'
            writer.writerow({'run_id':d['run_id'],'frame':frame.index,'time_ps':frame.time_fs/1000.,'observable':observable,'threshold':threshold,'fraction_basis':basis,**fractions,'mean_order':float(np.mean(val)),'max_order':float(np.max(np.abs(val)))})
    trajectory_qc={'status':'PASS','reasons':[],'warnings':[],'n_frames':nread,
        'max_force_eV_A':float(np.nanmax([x['fmax_eV_A'] for x in qc_rows])) if qc_rows else float('nan'),
        'minimum_AlN_A':float(np.nanmin([x.get('min_AlN_A',np.nan) for x in qc_rows])) if qc_rows else float('nan'),
        'minimum_ScN_A':float(np.nanmin([x.get('min_ScN_A',np.nan) for x in qc_rows])) if qc_rows else float('nan'),
        'max_BEC_abs':max((x['bec_abs_max'] for x in qc_rows if np.isfinite(x['bec_abs_max'])),default=float('nan')),
        'c_initial_A':qc_rows[0]['c_A'] if qc_rows else float('nan'),'c_final_A':qc_rows[-1]['c_A'] if qc_rows else float('nan')}
    for row in qc_rows:trajectory_qc['reasons'].extend(row['reasons']);trajectory_qc['warnings'].extend(row['warnings'])
    trajectory_qc['reasons']=sorted(set(trajectory_qc['reasons']));trajectory_qc['warnings']=sorted(set(trajectory_qc['warnings']))
    if trajectory_qc['reasons']:trajectory_qc['status']='FAILED'
    elif trajectory_qc['warnings']:trajectory_qc['status']='WARNING'
    overall='FAILED' if 'FAILED' in (thermo_qc['status'],trajectory_qc['status']) else 'WARNING' if 'WARNING' in (thermo_qc['status'],trajectory_qc['status']) else 'PASS'
    if overall=='FAILED' and not allow_failed:
        tmp.unlink(missing_ok=True);(out/'timeseries.csv').unlink(missing_ok=True)
    else:tmp.replace(out/'timeseries.csv')
    qc_doc={'run_id':d['run_id'],'status':overall,'thermo':thermo_qc,'trajectory':trajectory_qc}
    qc_fields=['source','frame','time_ps','temperature_K','potential_energy_eV','sigma_zz','c_A','max_force_eV_A','bec_abs_max','min_AlN_A','min_ScN_A','min_cation_cation_A','status','reason']
    with (out/'qc_timeseries.csv').open('w',newline='',encoding='utf-8') as qf:
        qw=csv.DictWriter(qf,fieldnames=qc_fields);qw.writeheader()
        if table is not None:
            values=table['values'];cols=table['columns'];interval=float(d.get('qc',{}).get('thermo_interval_fs',2.5))
            for i,row in enumerate(values):
                def tv(name):
                    ix=cols.get(name);return float(row[ix]) if ix is not None and ix<len(row) else float('nan')
                qw.writerow({'source':'thermo','frame':i,'time_ps':i*interval/1000.,'temperature_K':tv('temperature'),'potential_energy_eV':tv('potential_energy'),'sigma_zz':tv('sigma_zz'),'c_A':tv('c_A')})
        for row in qc_rows:
            reason='; '.join(row['reasons'])
            qw.writerow({'source':'trajectory','frame':row['frame'],'time_ps':row['time_ps'],'c_A':row['c_A'],'max_force_eV_A':row['fmax_eV_A'],'bec_abs_max':row['bec_abs_max'],
                         'min_AlN_A':row.get('min_AlN_A',float('nan')),'min_ScN_A':row.get('min_ScN_A',float('nan')),'min_cation_cation_A':row.get('min_cation_cation_A',float('nan')),'status':'FAILED' if reason else 'PASS','reason':reason})
    (out/'qc.json').write_text(json.dumps(qc_doc,indent=2,allow_nan=True),encoding='utf-8')
    source=Path(d['files']['trajectory']);elapsed=time.perf_counter()-started
    provenance={'analysis_version':__version__,'git_commit':os.environ.get('GIT_COMMIT','unknown'),'config_hash':cfg.digest,'source_run_id':d['run_id'],'source_file':str(source),'source_size_bytes':source.stat().st_size,'observable':observable,'threshold':threshold,'frame_stride':sampling,'frames_analyzed':nread,'elapsed_seconds':elapsed,'python':platform.python_version(),'platform':platform.platform()}
    (out/'analysis_manifest.json').write_text(json.dumps(provenance,indent=2),encoding='utf-8')
    if overall=='FAILED' and not allow_failed:raise RuntimeError(f"{d['run_id']} QC FAILED; downstream analysis skipped; see {out/'qc.json'}")
    return qc_doc

def analyze_project(project_path,outdir,backend='auto',allow_failed=False):
    project=load_yaml(project_path);base=Path(project_path).resolve().parent;out=Path(outdir);out.mkdir(parents=True,exist_ok=True);summary=[]
    for entry in project.get('runs',[]):
        mp=(base/entry).resolve() if not Path(entry).is_absolute() else Path(entry)
        run=load_run(mp);target=out/run.data['run_id']
        try:qc=analyze_run(mp,target,backend,allow_failed=allow_failed)
        except RuntimeError:
            qcp=target/'qc.json'
            if not qcp.exists():raise
            qc=json.loads(qcp.read_text())
        summary.append({'run_id':run.data['run_id'],'observable':run.data['analysis'].get('observable','q'),'status':qc['status'],'timeseries':str(target/'timeseries.csv') if (target/'timeseries.csv').exists() else ''})
    with (out/'project_summary.csv').open('w',newline='',encoding='utf-8') as f:
        writer=csv.DictWriter(f,fieldnames=['run_id','observable','status','timeseries']);writer.writeheader();writer.writerows(summary)
    return summary
