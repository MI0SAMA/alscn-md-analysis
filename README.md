# AlScN MD analysis (V0.1)

A memory-bounded analysis package for large GPUMD extended-XYZ trajectories. The C extension parses one frame at a time and returns compact typed buffers; Python/NumPy/ASE handle geometry, QC and outputs. Original simulation data is opened read-only.

## Install

```bash
python -m pip install -e '.[test]'
```

The editable install builds the native C reader. `backend=python` remains available for cross-checking.

## Run one trajectory

```bash
alscn-switch run examples/run_q.yaml --out analysis/A000_E0100
alscn-switch run examples/run_lambda.yaml --out analysis/A000_nucleus
```

`analysis.observable` is `q` or `lambda`. `q` dynamically rebuilds the four nearest N neighbors per frame and uses the tetrahedral c-axis third moment; default phases are q>0.2, q<−0.2 and transition otherwise. `lambda` projects minimum-image displacement onto matched parent/reversed reference endpoints and averages atomic values over the periodic 24×24 fractional a-b grid and classifies a cell at 0.5 (matching N0/N1). The lambda reference accepts the existing N0/N1 `reference.npz` keys `frac0` and `disp`, or a generic NPZ with `parent_positions`, `reversed_positions`, and `cell`; atom ordering must match the trajectory.

QC failures stop downstream analysis by default. `--allow-failed` preserves the computed timeseries for explicit diagnostics. Outputs are `qc.json`, `qc_timeseries.csv`, `timeseries.csv`, and `analysis_manifest.json`.

## Project matrix

```bash
alscn-switch project examples/project.yaml --out analysis/project
alscn-switch plot examples/project.yaml --analysis-dir analysis/project --out figures
alscn-switch benchmark /path/to/traj.xyz --frames 10
```

## Manifest and QC conventions

See `examples/run_q.yaml`. Relative trajectory/reference paths are resolved against the run YAML. GPUMD thermo columns default to temperature 0, kinetic energy 1, potential energy 2 and sigma_zz 8; c is measured from the trajectory cell. `thermo_interval_fs` defaults to 2.5 fs and can be overridden. Analysis provenance records source path/size, run ID, observable, threshold, config SHA256 and software version.

The native reader honors each frame's extxyz `Properties` declaration, including nonstandard attribute order. It parses species, position, unwrapped position, charge, force and BEC fields and rejects truncated/malformed frames. No trajectory is loaded in full.

## Scope

V0.1 includes trajectory/thermo QC, q/lambda observables, phase fractions, project summaries and matrix plots. Wall tracking, periodic cluster morphology, kinetic events, local chemistry, DFT export and full hysteresis analysis are later stages.


## Existing-data regression

The GPU source trajectories and reference assets stay outside this repository. To run the recorded leading-frame regressions, provide local paths to tiny readable trajectory samples (the first extxyz frame is sufficient) and the existing N1 `reference.npz`:

```bash
ALSCN_Q_TRAJECTORY=/path/to/A025_E0p050_first_frame.xyz \
ALSCN_N1_TRAJECTORY=/path/to/R14p4_w2p5_first_frame.xyz \
ALSCN_N1_REFERENCE=/path/to/reference.npz pytest -q tests/regression
```

The q reference checks the A025/E0.050 legacy dynamic-neighbor summary. The N1 reference checks the first `R14p4_w2p5` grid fraction (41 of 576 cells reversed). These tests are skipped unless paths are supplied.

## Direct source run on the GPU server

If the server Python lacks `venv/ensurepip`, the package can run from the checkout using its existing system dependencies:

```bash
cd /home/yao/otf_qnep/alscn-md-analysis
./scripts/build_native.sh
./scripts/alscn-switch run /path/to/run.yaml --out /path/to/separate-analysis-dir
```

The C extension is built for the active Python interpreter and stays ignored by Git. This path does not install packages into the server's system Python.
