#!/home/kwoksun2/anaconda3/bin/python3
import subprocess as sub
import numpy as np
import time
import os

NPOINT  = 8192 # 64**3
tmax    = 600.0

nT   = 16
nRho = 16


temp = np.logspace( 2, 3.5, nT )
rho  = np.logspace( 0, 14, nRho)
h2frac = np.logspace(-6, 0, 7 )

f_time = open("test_radau2a_time.txt", "w+")
f_input = open("test_radau2a_input.txt", "w+")
f_output = open("test_radau2a_output.txt", "w+")

COUNT = 0

for fH2 in h2frac:
  for T in temp:
    for Rho in rho:

      print("./test {0} {1} {2} {3} ".format( NPOINT, Rho, T, fH2 ))
      # output = sub.check_output( ["./test", str(NPOINT), str(Rho), str(T), str(fH2)  ], timeout=tmax , stderr=sub.STDOUT )
      try:
        tstart = time.time()
        output = sub.check_output( ["./test", str(NPOINT), str(Rho), str(T), str(fH2)  ], timeout=tmax )
        # output = sub.check_output( ["./test {0} {1} {2} {3} ".format( NPOINT, Rho, T, fH2 ) ], timeout=tmax )
        tend   = time.time()
        dt = tend - tstart
      except:
        dt = "NaN"
        # output = b"\n----------"
      COUNT += 1

      f_time.write("{0}\n".format(dt))
      f_input.write( ("{0} {1} {2}\n").format( Rho, T, fH2 ) )
      try:
        f_output.write("{0}\n".format(output.decode("utf-8")))
      except:
        f_output.write("Time: NaN\n")
      if (COUNT % 10 == 0):
        f_time.flush()
        f_input.flush()
        f_output.flush()
        os.fsync(f_time)
        os.fsync(f_input)
        os.fsync(f_output)

f_time.close()
f_input.close()
f_output.close()
