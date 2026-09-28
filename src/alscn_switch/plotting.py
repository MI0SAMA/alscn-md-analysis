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
        if not ts.exists():continue
        grouped.setdefault(cfg['analysis'].get('observable','q'),[]).append((cfg,ts))
    outputs=[]
    for observable,runs in grouped.items():
        n=len(runs);fig,axes=plt.subplots(n,1,figsize=(9,max(3.2*n,4)),squeeze=False,sharex=True)
        for ax,(cfg,path) in zip(axes[:,0],runs):
            with path.open(newline='') as f:rows=list(csv.DictReader(f))
            t=[float(x['time_ps']) for x in rows]
            for key,label,color in [('plus_fraction','+ phase','#c53929'),('minus_fraction','− phase','#315b9a'),('transition_fraction','transition','#888888')]:
                ax.plot(t,[float(x[key]) for x in rows],label=label,color=color,lw=1.25)
            st=cfg.get('structure',{});sim=cfg.get('simulation',{});name=st.get('name',cfg['run_id']);field=sim.get('field_VA')
            ax.set_title(f"{name}"+(f", E={field:g} V/Å" if isinstance(field,(int,float)) else ''),loc='left',fontsize=10)
            ax.set_ylabel('cation fraction');ax.set_ylim(-.03,1.03);ax.grid(alpha=.2);ax.legend(ncol=3,fontsize=8,loc='upper right')
        axes[-1,0].set_xlabel('time (ps)');fig.suptitle(f'AlScN phase fractions ({observable})');fig.tight_layout()
        dest=out/f'phase_matrix_{observable}.png';fig.savefig(dest,dpi=180);plt.close(fig);outputs.append(dest)
    return outputs
