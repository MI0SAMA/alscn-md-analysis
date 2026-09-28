"""Optional regressions against tiny leading-frame samples of GPU originals.

Set ALSCN_Q_TRAJECTORY, ALSCN_N1_TRAJECTORY and ALSCN_N1_REFERENCE to run;
these external data files are intentionally not stored in Git.
"""
import os
import numpy as np
import pytest
from alscn_switch.io.trajectory import iter_extxyz
from alscn_switch.polarization.order import compute_q,classify_q,phase_fractions,compute_lambda,classify_lambda_grid

@pytest.mark.skipif(not os.getenv('ALSCN_Q_TRAJECTORY'),reason='GPU trajectory path not configured')
def test_real_dynamic_q_matches_A025_E005_legacy_summary():
    frame=next(iter_extxyz(os.environ['ALSCN_Q_TRAJECTORY'],backend='c'));_,q=compute_q(frame)
    frac=phase_fractions(classify_q(q))
    assert abs(float(q.mean())-(-0.045361826))<2e-9
    assert abs(frac['plus_fraction']-0.4826049805)<1e-9
    assert abs(frac['minus_fraction']-0.5173950195)<1e-9

@pytest.mark.skipif(not (os.getenv('ALSCN_N1_TRAJECTORY') and os.getenv('ALSCN_N1_REFERENCE')),reason='GPU N1 paths not configured')
def test_real_N1_lambda_grid_matches_R14p4_w2p5_legacy_summary():
    frame=next(iter_extxyz(os.environ['ALSCN_N1_TRAJECTORY'],backend='c'))
    with np.load(os.environ['ALSCN_N1_REFERENCE'],allow_pickle=False) as z:reference={k:z[k] for k in z.files}
    lam=compute_lambda(frame,reference);phase,_,_=classify_lambda_grid(frame,lam)
    assert abs(frame.time_fs/1000-.02)<1e-10
    assert int(np.count_nonzero(phase==1))==41
    assert abs(float(np.mean(phase==1))-41/576)<1e-12
