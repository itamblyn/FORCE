#! /usr/bin/env python

import sys

if len(sys.argv) == 1:
     print 'usage: ' + sys.argv[0] + 'distance-force.int temperature'

import numpy

inputFile = open (sys.argv[1],'r')
INTEGRATED_array = []
#read line into array 
for line in inputFile.readlines():
    # add a new sublist
    INTEGRATED_array.append([])
    # loop over the elemets, split by whitespace
    for i in line.split():
        # convert to integer and append to the last
        # element of the list
        INTEGRATED_array[-1].append(float(i))
inputFile.close()

kB = 1.3806503E-23
T = float(sys.argv[2])

outputFile = open('RDF.' + str(int(T)) + '.dat','w')

for i in range(len(INTEGRATED_array)):
     outputFile.write(str(INTEGRATED_array[i][0]) + ' ' + str(numpy.exp(-INTEGRATED_array[i][1]/(kB*T))) + '\n')

outputFile.close()
