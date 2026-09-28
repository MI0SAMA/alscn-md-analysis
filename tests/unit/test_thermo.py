import numpy as np
from alscn_switch.qc.thermo import check_thermo,DEFAULT_COLUMNS

def test_runaway_fixture_reports_recorded_time_and_energy_jump():
    n=7000;a=np.zeros((n,19));a[:,0]=300.;a[:,2]=-250000.
    row=6567;a[row,0]=1114.;a[row,2]=-200000.
    result=check_thermo({'columns':DEFAULT_COLUMNS,'values':a},2.5,max_temperature_K=2000,natoms=32768)
    assert result['status']=='FAILED';assert abs(result['first_bad_time_ps']-16.4175)<1e-10
    assert 'temperature runaway (rolling median/MAD)' in result['reasons'];assert 'potential-energy discontinuity' in result['reasons']

def test_normal_table_passes():
    a=np.zeros((100,19));a[:,0]=300+np.sin(np.arange(100));a[:,2]=-250000-np.arange(100)
    result=check_thermo({'columns':DEFAULT_COLUMNS,'values':a})
    assert result['status']=='PASS'
