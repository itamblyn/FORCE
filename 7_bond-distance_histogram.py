#! /usr/local/bin/python

import sys

if len (sys.argv) == 1:
     print 'usage: ' + sys.argv[0] + ' distance.dat number_of_bond_length_bins [min_bond_length_value, max_bond_length_value]'

import numpy 


input_filename = sys.argv[1]
inputFile = open (input_filename,'r')
number_of_bond_length_bins = int(sys.argv[2])

infoline = inputFile.readline() # skip over info line

min_isolated_atom_distance = float(infoline.split()[5])
max_isolated_atom_distance = float(infoline.split()[7])

distance_array = []

for line in inputFile:

    distance_array.append(float(line.split()[0]))

inputFile.close()

if len(sys.argv) == 5:
     min_bond_length_value = float(sys.argv[3])
     max_bond_length_value = float(sys.argv[4])

else:
     min_bond_length_value = min(distance_array)
     max_bond_length_value = max(distance_array)

histogram = numpy.zeros(int(number_of_bond_length_bins + 1), dtype=numpy.float)

bond_length_bin_size = (max_bond_length_value - min_bond_length_value)/float(number_of_bond_length_bins)

for distance in distance_array:
     bin = int((distance - min_bond_length_value)/bond_length_bin_size)
     if bin > 0 and bin < len(histogram):
         histogram[bin] += 1

histogram /= numpy.sum(histogram)

output_filename = input_filename.rpartition('.')[0] + '.hist'
outputFile = open (output_filename, 'w')
outputFile.write('# distance occupancy, isolated atom distance( ' + str(min_isolated_atom_distance) + ' , ' + str(max_isolated_atom_distance) + ' )\n')

for bin_index in range(number_of_bond_length_bins):

     outputFile.write(repr(bin_index*bond_length_bin_size + bond_length_bin_size/2.0 + min_bond_length_value) + '  ')
     outputFile.write(str(histogram[bin_index]) + ' ')
     outputFile.write('\n')

outputFile.close()
