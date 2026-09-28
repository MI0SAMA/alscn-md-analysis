import numpy as np
from alscn_switch.io.trajectory import Frame
from alscn_switch.structure.neighbors import dynamic_cation_n_neighbors,minimum_distances

def test_neighbor_search_crosses_periodic_boundary_in_sheared_cell():
    cell=np.array([[10.,0,0],[1.,10.,0],[0,0,10.]])
    # Four N sites surround a Sc; one nearest neighbor is across x PBC.
    pos=np.array([[.2,5.,5.],[9.4,5.,4.8],[.8,5.,5.2],[.2,4.4,5.2],[.2,5.6,4.8]])
    fr=Frame(1,0,np.array([2,3,3,3,3],np.uint8),pos,pos.copy(),cell,np.full(5,np.nan),np.full((5,3),np.nan),np.full((5,9),np.nan))
    c,n,d=dynamic_cation_n_neighbors(fr)
    assert c.shape==(1,4);assert set(n[0])=={1,2,3,4};assert np.max(np.linalg.norm(d[0],axis=1))<1.5
    assert minimum_distances(fr)['min_ScN_A']<1.
