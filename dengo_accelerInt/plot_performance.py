#!/home/kwoksun2/anaconda3/bin/python
import os
import numpy
import matplotlib.pyplot as plt
import subprocess
import numpy as np

tout = subprocess.check_output(["grep", "-r","Time:", "test_radau2a_output.txt"])

t = [ numpy.float(str(i)) for i in tout.decode("utf-8").replace("Time:","").split("sec\n") if i != ""]

t_klu = []
with open( "/home/kwoksun2/dengo-merge/dengo_data_structure/test_klu_time.txt" ) as f:
  for l in f.readlines():
    t_klu.append(float(l))

t = numpy.array(t)
t_klu = numpy.array(t_klu)

lbin = numpy.logspace(-3, 3, 300)
print( numpy.isnan(t).sum() )
plt.hist(t[~numpy.isnan(t)], bins = lbin, label='gpu relto=1e-6')
plt.hist(t_klu, bins=lbin, alpha=0.5, label='cpu klu solver')
plt.xscale("log")
plt.ylabel('count')
plt.legend()
plt.xlabel("runtime (s) / 8192 cells")
plt.savefig("runtime_count.png")


lbin = numpy.logspace(-2.5, 2.5, 100)
print( numpy.isnan(t).sum() )
t_ratio = t_klu / t
plt.clf()
plt.hist( t_ratio[~numpy.isnan(t_ratio)], bins=lbin, label='t_klu / t' )
plt.xscale("log")
plt.ylabel('count')
plt.legend()
plt.xlabel("t_klu/t")
plt.savefig("runtime_comparison_count.png")



nT   = 16
nRho = 16


temp = np.logspace( 2, 3.5, nT )
rho  = np.logspace( 0, 14, nRho)
h2frac = np.logspace(-6, 0, 7 )

i = 6

dtemp = np.log10(temp[1]/ temp[0])/2.0
drho  = np.log10(rho[1] / rho[0])/2.0

X = np.logspace( 0-drho, 14+drho, nRho+1 )
y = np.logspace( 2-dtemp, 3.5+dtemp, nT + 1)

t_reshape = t.reshape( 7, 16, 16)
t_klu_reshape = t_klu.reshape( 7,16,16  )
for fh2 in h2frac:
  _t =     t[i*256: (i+1)*256].reshape(16,16)
  _t_klu = t_klu[i*256: (i+1)*256].reshape(16,16)

  print(  t_reshape[i] - _t)
  print( t_klu_reshape[i] - _t_klu )

  plt.clf()
  im = plt.pcolormesh(X, y, np.log10(_t_klu / _t) , \
             vmin=-2, vmax = 2, cmap = "coolwarm_r")

  cbar = plt.colorbar( im , fraction=0.046, pad = 0.04)
  cbar.ax.set_ylabel("log ($t_{cpu}$ / $t_{gpu}$)")
  plt.xscale("log")
  plt.yscale("log")
  plt.xlabel("rho (cm**-3)")
  plt.ylabel("temperature (K)")
  plt.savefig("fh_1e-{}_2d.png".format(i))
  i -= 1

import h5py
f = h5py.File("dengo_performance.h5", 'w')
f.create_dataset( "H2_fraction_0", data = h2frac )
f.create_dataset( "rho_2", data = X )
f.create_dataset( "temp_1", data = y  )
f.create_dataset( "t_cpu", data=t_klu_reshape)
f.create_dataset( "t_gpu" , data=t_reshape )
f.close()

