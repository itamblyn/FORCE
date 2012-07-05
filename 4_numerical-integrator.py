#! /usr/bin/env python

import sys

if len(sys.argv) == 1:
     print 'usage: ' + sys.argv[0] + ' distance-force.hist'

inputFile = open (sys.argv[1],'r')
inputFile.readline()
AVERAGED_array = []
#read line into array 
for line in inputFile.readlines():
    # add a new sublist
    AVERAGED_array.append([])
    # loop over the elemets, split by whitespace
    for i in line.split():
        # convert to integer and append to the last
        # element of the list
        AVERAGED_array[-1].append(float(i))
inputFile.close()

bin_size = (AVERAGED_array[1][0] - AVERAGED_array[0][0])*0.529177

outputFile = open('distance-force.int', 'w')

for i in range(len(AVERAGED_array)):
     value = 0
     for j in range(i):
          value += AVERAGED_array[j][1]*bin_size*4.35974417E-18
     outputFile.write(str(AVERAGED_array[i][0]) + ' ' + str(-value) + '\n')

outputFile.close()
