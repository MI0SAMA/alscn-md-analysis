"""Project-level time-series and phase-fraction matrix plots."""
from __future__ import annotations
import csv
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import yaml

def plot_project(project_path,analysis_dir,outdir):
    project_path=Path(project_path).resolve();base=project_path.parent;data=yaml.safe_load(project_path.read_text());root=Path(analysis_dir);out=Path(outdir);out.mkdir(parents=True,exist_ok=True)
    grouped={}
    for manifest in data.get('runs',[]):
        mp=(base/manifest).resolve();cfg=yaml.safe_load(mp.read_text());rid=cfg['run_id'];ts=root/rid/'timeseries.csv'
        if ts.exists():grouped.setdefault(cfg['analysis'].get('observable','q'),[]).append((cfg,ts))
    outputs=[]
    for observable,runs in grouped.items():
        motifs=list(dict.fromkeys(str(c.get('structure',{}).get('name',c['run_id'])) for c,_ in runs))
        fields=sorted({c.get('simulation',{}).get('field_VA') for c,_ in runs if isinstance(c.get('simulation',{}).get('field_VA'),(int,float))})
        if not fields:fields=[None]
        fig,axes=plt.subplots(len(motifs),len(fields),figsize=(4.1*len(fields),3.0*len(motifs)),squeeze=False,sharex=True,sharey=True)
        for cfg,path in runs:
            name=str(cfg.get('structure',{}).get('name',cfg['run_id']));field=cfg.get('simulation',{}).get('field_VA')
            i=motifs.index(name);j=fields.index(field) if field in fields else 0;ax=axes[i,j]
            with path.open(newline='') as f:rows=list(csv.DictReader(f))
            t=[float(x['time_ps']) for x in rows]
            for key,label,color in [('plus_fraction','+ phase','#c53929'),('minus_fraction','− phase','#315b9a'),('transition_fraction','transition','#888888')]:
                ax.plot(t,[float(x[key]) for x in rows],label=label,color=color,lw=1.25)
            ax.grid(alpha=.2);ax.set_ylim(-.03,1.03)
        for i,name in enumerate(motifs):axes[i,0].set_ylabel(f'{name}\nfraction')
        for j,field in enumerate(fields):axes[0,j].set_title('all fields' if field is None else f'E = {field:g} V/Å')
        for ax in axes[-1,:]:ax.set_xlabel('time (ps)')
        handles,labels=axes[0,0].get_legend_handles_labels()
        if handles:fig.legend(handles,labels,loc='upper center',ncol=3,bbox_to_anchor=(.5,1.01),fontsize=9)
        fig.suptitle(f'AlScN phase fractions ({observable})');fig.tight_layout(rect=(0,0,1,.97))
        dest=out/f'phase_matrix_{observable}.png';fig.savefig(dest,dpi=180);plt.close(fig);outputs.append(dest)
    return outputs
