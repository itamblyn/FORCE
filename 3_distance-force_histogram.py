#! /usr/local/bin/python

import sys

if len(sys.argv) == 1:
     print 'usage: ' + sys.argv[0] + ' number_bond_length_bins number_of_atom_distance_bins'

import numpy

inputFile = open ('distance-force.dat','r')
number_of_bond_length_bins = int(sys.argv[1])
number_of_atom_distance_bins = int(sys.argv[2])

min_bond_length_value = 1E6
max_bond_length_value = 0
min_atom_distance_value = 1E6
max_atom_distance_value = 0

inputFile.readline() # skip over info line

input_array = []

for line in inputFile:
    input_array.append([])
    for i in line.split():
        input_array[-1].append(float(i))
    min_bond_length_value   = min(input_array[-1][0], min_bond_length_value)
    max_bond_length_value   = max(input_array[-1][0], max_bond_length_value)
    min_atom_distance_value = min(input_array[-1][2], min_atom_distance_value)
    max_atom_distance_value = max(input_array[-1][2], max_atom_distance_value)

inputFile.close()

# print 'input read complete'

bond_length_bin_size   = (max_bond_length_value   - min_bond_length_value)  /float(number_of_bond_length_bins)
atom_distance_bin_size = (max_atom_distance_value - min_atom_distance_value)/float(number_of_atom_distance_bins)

if atom_distance_bin_size == 0:
    atom_distance_bin_size = 1

force_array = []
distance_array = []

for i in range(number_of_atom_distance_bins + 1):
    force_array.append([])
    distance_array.append([])
    for j in range(number_of_bond_length_bins + 1):    # this creates empty force array
        force_array[i].append([])

for element in input_array:

    bond_length   = element[0]
    force         = element[1]
    atom_distance = element[2]

    atom_distance_bin        = int((atom_distance - min_atom_distance_value)/atom_distance_bin_size)
    bond_length_distance_bin = int((bond_length -   min_bond_length_value)/bond_length_bin_size)

    force_array[atom_distance_bin][bond_length_distance_bin].append(force)
    distance_array[atom_distance_bin].append(bond_length)


for i in range(number_of_atom_distance_bins):
    outputFile = open('distance.' + str(i) + '.dat', 'w')
    outputFile.write('# distance, atom distance ( ' + str(min_atom_distance_value + i*atom_distance_bin_size) + ' , ' + str(min_atom_distance_value + (i + 1)*atom_distance_bin_size) + ' )\n')
    for distance in distance_array[i]:
        outputFile.write(str(distance) + '\n')
    outputFile.close()

for i in range(number_of_atom_distance_bins):

    outputFile = open('distance-force.' + str(i) + '.hist', 'w')
    outputFile.write('# distance mean std_err, isolated atom distance ( ' + str(min_atom_distance_value + i*atom_distance_bin_size) + ' , ' + str(min_atom_distance_value + (i + 1)*atom_distance_bin_size) + ' )\n')

    for j in range(len(force_array[i])):
        if (len(force_array[i][j]) > 0):

             bond_length_bin_value = j*bond_length_bin_size + min_bond_length_value + bond_length_bin_size/2.0
             outputFile.write(repr(bond_length_bin_value))

             force_array[i][j] = numpy.array(force_array[i][j])   # converts each sublist to a numpy array

             mean = numpy.mean(force_array[i][j])
             standard_deviation = numpy.std(force_array[i][j])
             standard_error = standard_deviation/numpy.sqrt(len(force_array[i][j]))

             outputFile.write(' ' + repr(mean))
             outputFile.write(' ' + repr(standard_error))
             outputFile.write('\n')

    outputFile.close()
