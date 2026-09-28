import numpy as np
import pytest
from alscn_switch.io.trajectory import iter_extxyz

def write_xyz(path,properties='species:S:1:charge:R:1:pos:R:3:unwrapped_position:R:3:forces:R:3:bec:R:9'):
    rows=[]
    # Deliberately non-default column order; includes all native properties.
    for el,pos,q in [('Sc',(0.2,5,5),1.1),('N',(9.4,5,4.8),-1.0),('N',(0.8,5,5.2),-1.0),('N',(0.2,4.4,5.2),-1.0),('N',(0.2,5.6,4.8),-1.0)]:
        fields=[el,str(q),*(str(x) for x in pos),*(str(x) for x in pos), '0','0','0',*(['1']*9)]
        rows.append(' '.join(fields))
    path.write_text(f'5\nTime=10.0 Lattice="10 0 0 0 10 0 0 0 10" Properties={properties} pbc="T T T"\n'+'\n'.join(rows)+'\n')

def test_c_and_python_reader_handle_properties_order(tmp_path):
    p=tmp_path/'small.xyz';write_xyz(p)
    c=list(iter_extxyz(p,backend='c'))[0];py=list(iter_extxyz(p,backend='python'))[0];ase=list(iter_extxyz(p,backend='ase'))[0]
    np.testing.assert_array_equal(c.species,py.species);np.testing.assert_allclose(c.positions,py.positions)
    np.testing.assert_allclose(c.charges,py.charges);np.testing.assert_allclose(c.forces,py.forces);np.testing.assert_allclose(c.positions,ase.positions);np.testing.assert_array_equal(c.species,ase.species);np.testing.assert_allclose(c.unwrapped_positions,ase.unwrapped_positions);np.testing.assert_allclose(c.bec,ase.bec);np.testing.assert_allclose(c.forces,ase.forces);assert c.time_fs==10

def test_truncated_frame_is_explicit_error(tmp_path):
    p=tmp_path/'bad.xyz';p.write_text('2\nTime=0 Lattice="10 0 0 0 10 0 0 0 10" Properties=species:S:1:pos:R:3\nAl 0 0 0\n')
    with pytest.raises(ValueError,match='truncated'):
        list(iter_extxyz(p,backend='c'))
