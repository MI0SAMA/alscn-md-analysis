"""Command-line interface."""
from __future__ import annotations
import argparse,json
from pathlib import Path
from .analysis import analyze_run,analyze_project
from .io.trajectory import benchmark_readers
from .plotting import plot_project

def main(argv=None):
    parser=argparse.ArgumentParser(prog='alscn-switch',description='Streaming AlScN MD trajectory QC and phase analysis')
    sub=parser.add_subparsers(dest='command',required=True)
    run=sub.add_parser('run',help='analyze one run manifest');run.add_argument('manifest');run.add_argument('--out',required=True);run.add_argument('--backend',choices=['auto','c','python'],default='auto');run.add_argument('--stride',type=int);run.add_argument('--allow-failed',action='store_true')
    project=sub.add_parser('project',help='analyze all runs in project YAML');project.add_argument('manifest');project.add_argument('--out',required=True);project.add_argument('--backend',choices=['auto','c','python'],default='auto');project.add_argument('--allow-failed',action='store_true')
    plot=sub.add_parser('plot',help='plot analyzed project time series');plot.add_argument('manifest');plot.add_argument('--analysis-dir',required=True);plot.add_argument('--out',required=True)
    bench=sub.add_parser('benchmark',help='compare C and Python extxyz readers');bench.add_argument('trajectory');bench.add_argument('--frames',type=int,default=10)
    args=parser.parse_args(argv)
    if args.command=='run':
        try:result=analyze_run(args.manifest,args.out,args.backend,args.stride,args.allow_failed)
        except RuntimeError as e:print(e);return 2
        print(json.dumps(result,indent=2));return 0
    if args.command=='project':
        rows=analyze_project(args.manifest,args.out,args.backend,args.allow_failed);print(json.dumps(rows,indent=2));return 0
    if args.command=='plot':
        files=plot_project(args.manifest,args.analysis_dir,args.out);print('\n'.join(map(str,files)));return 0
    if args.command=='benchmark':
        print(json.dumps(benchmark_readers(args.trajectory,args.frames),indent=2));return 0
    return 2
if __name__=='__main__':raise SystemExit(main())
