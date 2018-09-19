#/home/kwoksun2/anaconda2/bin/python

from argparse import ArgumentParser
import h5py
import numpy

parser = ArgumentParser( description = "Compute the relative error of the input file to the standard solution file" )
parser.add_argument("-t", "--test-file", dest="test_file", help = " input file to be tested agaist ")
parser.add_argument("-s", "--solution-file", dest="solution_file", help = " standard solution with higher accuracy threshold ")

args = parser.parse_args()

tf = h5py.File(args.test_file)
sf = h5py.File(args.solution_file)

error = []
for k in tf.keys():
    error = ( numpy.abs(tf[k].value[0] - sf[k].value[0]) / (sf[k].value[0]) )
    print("{}_error: {}; {}; {}".format(k ,error,tf[k].value[0], sf[k].value[0]  ) )

tf.close()
sf.close()


#for _ in range(args.n):
#    print("Hello, {}".format(args.user_name))
