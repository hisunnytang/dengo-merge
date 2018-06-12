import numpy
import h5py

f = h5py.File('oldIC.h5')
data = {}

size = 128

for key in f.keys():
    if key not in ['T', 'ge']:
        data[key] = f[key]*numpy.ones((size))
    else:
        data[key] = f[key]*numpy.ones((size))
    print(key,f[key].value)


f.close()

f = h5py.File('eq_H2_2_initial_conditions.h5')
for k,v in data.items():
    f.create_dataset(k, data= v)
f.close()
