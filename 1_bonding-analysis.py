#! /usr/bin/env python

import sys

if len(sys.argv) == 1:
    print 'usage: ' + sys.argv[0] + ' FTRAJECTORY cell_size natoms position_force_bonding_file'

import numpy
from progressBar import *

position_force_file = sys.argv[1]
cell_size = float(sys.argv[2])
number_of_particles = int(sys.argv[3])
position_force_bonding_file = sys.argv[4]


def pbc_round(input_value):
     i = int(input_value)
     if (abs(input_value-i) >= 0.5):
          if (input_value > 0): i+= 1
          if (input_value < 0): i-= 1
     return i

inputFile = open (position_force_file,'r')
POSITION_FORCE_array = []
#read line into array 
for line in inputFile.readlines():
 ##########################
##                        ##
##  xx yy zz fx fy fz b?  ##
##                        ##
 ##########################
    # add a new sublist
    POSITION_FORCE_array.append([])
    POSITION_FORCE_array[-1].append(float(line.split()[1]))
    POSITION_FORCE_array[-1].append(float(line.split()[2]))
    POSITION_FORCE_array[-1].append(float(line.split()[3]))
    POSITION_FORCE_array[-1].append(float(line.split()[7]))
    POSITION_FORCE_array[-1].append(float(line.split()[8]))
    POSITION_FORCE_array[-1].append(float(line.split()[9]))
    POSITION_FORCE_array[-1].append(-10)
inputFile.close()

number_of_snapshots = int(len(POSITION_FORCE_array)/float(number_of_particles))

CELL_array = cell_size*numpy.ones((number_of_snapshots,3), dtype=numpy.float)

print 'input read complete'

sanity = (CELL_array[0][0]**2 + CELL_array[0][1]**2 + CELL_array[0][2]**2)**(1.0/2.0)

progess = progressBar(0, number_of_snapshots, 77)

s = 0 # counts over snapshots

while s < number_of_snapshots:

     molecule_number = 0
   
     p = 0 # counts over particles
     
     while p < number_of_particles:
          
          if (POSITION_FORCE_array[s*number_of_particles + p][6] == -10):
     
               o = 0 # counts over other particles

               minimum_distance = sanity
          
               closest_to_p = p

               while o < number_of_particles:

                    dx = POSITION_FORCE_array[s*number_of_particles + p][0] - POSITION_FORCE_array[s*number_of_particles + o][0]
                    dy = POSITION_FORCE_array[s*number_of_particles + p][1] - POSITION_FORCE_array[s*number_of_particles + o][1]
                    dz = POSITION_FORCE_array[s*number_of_particles + p][2] - POSITION_FORCE_array[s*number_of_particles + o][2]

                    dx -= CELL_array[s][0]*pbc_round(dx/CELL_array[s][0])
                    dy -= CELL_array[s][1]*pbc_round(dy/CELL_array[s][1])
                    dz -= CELL_array[s][2]*pbc_round(dz/CELL_array[s][2])

                    distance = (dx**2 + dy**2 + dz**2)**(0.5)

                    if (distance > sanity): print 'Warning, problem with pbc'
               
                    if (distance < minimum_distance and distance != 0.0):
                         minimum_distance = distance
                         closest_to_p = o 
           
                    o +=1

               minimum_distance = sanity

               o = closest_to_p
               oo = 0

               if (POSITION_FORCE_array[s*number_of_particles + o][6] == -10):
          
                    while oo < number_of_particles:
               
                         dx = POSITION_FORCE_array[s*number_of_particles + o][0] - POSITION_FORCE_array[s*number_of_particles + oo][0]
                         dy = POSITION_FORCE_array[s*number_of_particles + o][1] - POSITION_FORCE_array[s*number_of_particles + oo][1]
                         dz = POSITION_FORCE_array[s*number_of_particles + o][2] - POSITION_FORCE_array[s*number_of_particles + oo][2]

                         dx -= CELL_array[s][0]*pbc_round(dx/CELL_array[s][0])
                         dy -= CELL_array[s][1]*pbc_round(dy/CELL_array[s][1])
                         dz -= CELL_array[s][2]*pbc_round(dz/CELL_array[s][2])

                         distance = (dx**2 + dy**2 + dz**2)**(0.5)
                         
                         if (distance > sanity): print 'Warning, problem with pbc'
               
                         if (distance < minimum_distance and distance != 0.0):
                              minimum_distance = distance
                              closest_to_o = oo 
          
                         oo +=1

                    if (closest_to_p == o and closest_to_o == p):
                         molecule_number += 1
                         POSITION_FORCE_array[s*number_of_particles + p][6] = o
                         POSITION_FORCE_array[s*number_of_particles + o][6] = p
                         
                    else: POSITION_FORCE_array[s*number_of_particles + p][6] = -1
                    
               else: POSITION_FORCE_array[s*number_of_particles + p][6] = -1                    

          p +=1
     
     progess.updateAmount(s)
     print progess, "\r",

     s += 1


 ####################################
##                                  ##
## bonding analysis is now complete ##
##                                  ##
 ####################################

position_force_bonding = open(position_force_bonding_file, 'w')

for atoms in POSITION_FORCE_array:
     for element in atoms:
          position_force_bonding.write(str(element) + ' ')
     position_force_bonding.write('\n')

position_force_bonding.close()

print
