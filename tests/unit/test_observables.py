import numpy as np
from alscn_switch.polarization.order import classify_q,classify_lambda,phase_fractions,compute_lambda
from alscn_switch.io.trajectory import Frame

def test_classifiers_and_fractions():
    np.testing.assert_array_equal(classify_q(np.array([.21,-.21,.2,-.2,0.])),[1,-1,0,0,0])
    np.testing.assert_array_equal(classify_lambda(np.array([.49,.5,.51])),[-1,-1,1])
    assert phase_fractions(np.array([1,-1,0]))=={'plus_fraction':1/3,'minus_fraction':1/3,'transition_fraction':1/3}

def test_lambda_uses_minimum_image_fractional_displacement():
    cell=np.diag([10.,10.,10.]);parent=np.array([[9.8,5.,5.]]);rev=np.array([[.8,5.,5.]])
    # Current at x=0.3 is halfway along the short periodic path 9.8 -> 0.8.
    frame=Frame(1,0.,np.array([1],np.uint8),np.array([[.3,5.,5.]]),np.array([[.3,5.,5.]]),cell,np.array([np.nan]),np.full((1,3),np.nan),np.full((1,9),np.nan))
    val=compute_lambda(frame,{'parent_positions':parent,'reversed_positions':rev,'cell':cell})
    np.testing.assert_allclose(val[0],.5,atol=1e-12)

def test_lambda_grid_uses_cell_average_not_atom_count():
    from alscn_switch.polarization.order import classify_lambda_grid
    cell=np.diag([10.,10.,10.]);pos=np.array([[1.,1.,1.],[1.1,1.1,1.],[6.,6.,1.]])
    frame=Frame(1,0.,np.array([1,2,3],np.uint8),pos,pos.copy(),cell,np.full(3,np.nan),np.full((3,3),np.nan),np.full((3,9),np.nan))
    phase,mean,count=classify_lambda_grid(frame,np.array([.2,.8,.9]),.5,(2,2))
    assert phase[0,0]==-1 and count[0,0]==2
    assert phase[1,1]==1 and count[1,1]==1
