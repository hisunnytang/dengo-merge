import numpy
import h5py

f = h5py.File('openmp_example/cvdls_9species_initial_conditions.h5')
data = {}

size = 1024

for key in f.keys():
    if key not in ['T', 'ge']:
        data[key] = f[key]*numpy.logspace(0,0.1,size)
    else:
        data[key] = f[key]*numpy.ones((size))

f.close()

f = h5py.File('newIC.h5')
for k,v in data.items():
    f.create_dataset(k, data= v)
f.close()
