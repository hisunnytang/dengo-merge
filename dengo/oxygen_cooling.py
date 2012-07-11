import h5py
import numpy as na
from primordial_cooling import CoolingRate, CoolingAction, cooling_rates_table, cooling_action_table

CoolingRate.init_temperature((1.0, 1e9))
T = CoolingRate.T
logT = CoolingRate.logT
tev= CoolingRate.tev
logtev = CoolingRate.logtev

# oxygen cooling rates
f = h5py.File('o_ion_by_ion_cooling.h5')
data = f['Table']
for i in range(len(data.dtype.names) - 2):
    vals = na.interp(T, data['T'], data['o_%i' %(i+1)])
    cooling_rates_table['o_%i_c' %(i+1)] = CoolingRate('o_%i_c' %(i+1), vals)
    cooling_action_table["o_%i_c" %(i+1)] = CoolingAction(
        ["o_%i_c" %(i+1)], ["o_%i" %(i+1),"de"], "-o_%i_c * o_%i * de" %(i+1,i+1))
