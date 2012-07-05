#! /usr/local/bin/python

import numpy, matplotlib, pylab, sys

if len(sys.argv) == 1:
     print 'usage: ' + sys.argv[0] + ' distance.dat number_of_bins'


input_filename = sys.argv[1]

inputFile = open (input_filename,'r')

number_of_bins = int(sys.argv[2])

inputFile.readline() # skip over info line
DISTANCE_array = []
for line in inputFile:
     DISTANCE_array.append(float(line.split()[0]))
inputFile.close()

png = input_filename.rpartition('.')[0] + '.png'


DISTANCE_array = numpy.reshape(DISTANCE_array,( len(DISTANCE_array) ))

pylab.hist(DISTANCE_array, number_of_bins, 1)
pylab.axis([1.0,2.5,0,5])

pylab.savefig(png)
