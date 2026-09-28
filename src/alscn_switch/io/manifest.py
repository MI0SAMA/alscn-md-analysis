"""YAML manifest loading and reproducibility metadata."""
from __future__ import annotations
from dataclasses import dataclass
import hashlib, json
from pathlib import Path
from typing import Any
import yaml

@dataclass(slots=True)
class RunConfig:
    path: Path
    data: dict[str, Any]
    digest: str

def load_yaml(path: str | Path) -> dict[str, Any]:
    p=Path(path).resolve(); data=yaml.safe_load(p.read_text(encoding='utf-8')) or {}
    if not isinstance(data,dict): raise ValueError(f'expected mapping at {p}')
    return data

def load_run(path: str | Path) -> RunConfig:
    p=Path(path).resolve(); data=load_yaml(p)
    missing=[k for k in ('run_id','files','analysis') if k not in data]
    if missing: raise ValueError(f'run manifest missing: {", ".join(missing)}')
    if 'trajectory' not in data['files']: raise ValueError('files.trajectory is required')
    for key,val in list(data['files'].items()):
        if isinstance(val,str) and not Path(val).is_absolute(): data['files'][key]=str((p.parent/val).resolve())
    if isinstance(data['analysis'].get('reference'),str) and not Path(data['analysis']['reference']).is_absolute():
        data['analysis']['reference']=str((p.parent/data['analysis']['reference']).resolve())
    canonical=json.dumps(data,sort_keys=True,separators=(',',':'),default=str)
    return RunConfig(p,data,hashlib.sha256(canonical.encode()).hexdigest())
