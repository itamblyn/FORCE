#! /usr/local/bin/python

import sys

if len(sys.argv) == 1:
     print 'usage: ' + sys.argv[0] + ' position_force_bonding_file cell_size natoms'

import numpy
from progressBar import *

position_force_bonding_file = sys.argv[1]
cell_size = float(sys.argv[2])
number_of_particles = int(sys.argv[3])

def pbc_round(input_value):
     i = int(input_value)
     if (abs(input_value-i) >= 0.5):
          if (input_value > 0): i+= 1
          if (input_value < 0): i-= 1
     return i

inputFile = open (position_force_bonding_file,'r')
POSITION_FORCE_array = []
#read line into array 
for line in inputFile:
    POSITION_FORCE_array.append([])
    POSITION_FORCE_array[-1].append(float(line.split()[0])) # xx
    POSITION_FORCE_array[-1].append(float(line.split()[1])) # yy
    POSITION_FORCE_array[-1].append(float(line.split()[2])) # zz
    POSITION_FORCE_array[-1].append(float(line.split()[3])) # fx
    POSITION_FORCE_array[-1].append(float(line.split()[4])) # fy
    POSITION_FORCE_array[-1].append(float(line.split()[5])) # fz
    POSITION_FORCE_array[-1].append(int(line.split()[6]))   # bonding?
inputFile.close()

number_of_snapshots = int(len(POSITION_FORCE_array)/float(number_of_particles))

CELL_array = cell_size*numpy.ones((number_of_snapshots,3), dtype=numpy.float)

print 'input read complete'

outputFile = open('distance-force.dat', 'w')
outputFile.write('# bond_length force distance_to_nearest_atom\n')

progess = progressBar(0, number_of_snapshots, 77)

for s in range(number_of_snapshots):

     atom_list = []
     
     for p in range(number_of_particles):
     
          bonding = POSITION_FORCE_array[s*number_of_particles + p][6]

          if (bonding == -1):
               atom_list.append(p)

     for p in range(number_of_particles):

          bonding = POSITION_FORCE_array[s*number_of_particles + p][6]

          if (bonding != -1):
                
               px = POSITION_FORCE_array[s*number_of_particles + p][0]            # these are the coordinates of the particle
               py = POSITION_FORCE_array[s*number_of_particles + p][1]
               pz = POSITION_FORCE_array[s*number_of_particles + p][2]

               o  = POSITION_FORCE_array[s*number_of_particles + p][6]            # this is the index of the particle to which it is bonded

               for oo in range(number_of_particles):                                     # now shift the coordinates of the array so that the original atom is at the origin, and pbc
                    x = POSITION_FORCE_array[s*number_of_particles + oo][0] - px
                    y = POSITION_FORCE_array[s*number_of_particles + oo][1] - py
                    z = POSITION_FORCE_array[s*number_of_particles + oo][2] - pz
                    POSITION_FORCE_array[s*number_of_particles + oo][0] = x -(int(x/CELL_array[s][0]+number_of_snapshots+0.5)-number_of_snapshots)*CELL_array[s][0]
                    POSITION_FORCE_array[s*number_of_particles + oo][1] = y -(int(y/CELL_array[s][1]+number_of_snapshots+0.5)-number_of_snapshots)*CELL_array[s][1]
                    POSITION_FORCE_array[s*number_of_particles + oo][2] = z -(int(z/CELL_array[s][2]+number_of_snapshots+0.5)-number_of_snapshots)*CELL_array[s][2]

               ox = POSITION_FORCE_array[s*number_of_particles + o][0]            # these are now the coordinates of the atom to which it is bonded
               oy = POSITION_FORCE_array[s*number_of_particles + o][1]
               oz = POSITION_FORCE_array[s*number_of_particles + o][2]

               for oo in range(number_of_particles):              # put the MOLECULE at the centre of the cell, and reapply pbc

                    x = POSITION_FORCE_array[s*number_of_particles + oo][0] - ox/2.0
                    y = POSITION_FORCE_array[s*number_of_particles + oo][1] - oy/2.0
                    z = POSITION_FORCE_array[s*number_of_particles + oo][2] - oz/2.0
                    POSITION_FORCE_array[s*number_of_particles + oo][0] = x -(int(x/CELL_array[s][0]+number_of_snapshots+0.5)-number_of_snapshots)*CELL_array[s][0]
                    POSITION_FORCE_array[s*number_of_particles + oo][1] = y -(int(y/CELL_array[s][1]+number_of_snapshots+0.5)-number_of_snapshots)*CELL_array[s][1]
                    POSITION_FORCE_array[s*number_of_particles + oo][2] = z -(int(z/CELL_array[s][2]+number_of_snapshots+0.5)-number_of_snapshots)*CELL_array[s][2]

               px = POSITION_FORCE_array[s*number_of_particles + p][0] # these are now the shifted/pbc coordinates
               py = POSITION_FORCE_array[s*number_of_particles + p][1]
               pz = POSITION_FORCE_array[s*number_of_particles + p][2]
               ox = POSITION_FORCE_array[s*number_of_particles + o][0]
               oy = POSITION_FORCE_array[s*number_of_particles + o][1]
               oz = POSITION_FORCE_array[s*number_of_particles + o][2]

               bond_length = ((px - ox)**2 + (py - oy)**2 + (pz - oz)**2)**(1.0/2.0)

               bx = POSITION_FORCE_array[s*number_of_particles + o][0] - POSITION_FORCE_array[s*number_of_particles + p][0]
               by = POSITION_FORCE_array[s*number_of_particles + o][1] - POSITION_FORCE_array[s*number_of_particles + p][1]
               bz = POSITION_FORCE_array[s*number_of_particles + o][2] - POSITION_FORCE_array[s*number_of_particles + p][2]

               fx = POSITION_FORCE_array[s*number_of_particles + p][3]
               fy = POSITION_FORCE_array[s*number_of_particles + p][4]
               fz = POSITION_FORCE_array[s*number_of_particles + p][5]

               force = -(1.0/(bond_length))*(fx*bx + fy*by + fz*bz)
               
               atom_distance_list = []

               for atom in atom_list:
                    
                    ax = POSITION_FORCE_array[s*number_of_particles + atom][0]
                    ay = POSITION_FORCE_array[s*number_of_particles + atom][1]
                    az = POSITION_FORCE_array[s*number_of_particles + atom][2]
          
                    atom_distance = (ax**2 + ay**2 + az**2)**(1.0/2.0)
                    atom_distance_list.append(atom_distance)
               
               if len(atom_distance_list) > 0:
                    distance_to_nearest_atom = min(atom_distance_list)
               else:
                    distance_to_nearest_atom = ((CELL_array[s][0]/2.0)**2 + (CELL_array[s][1]/2.0)**2 + (CELL_array[s][2]/2.0)**2)**(1.0/2.0)

               outputFile.write(str(bond_length) + ' ' + str(force) + ' ' + str(distance_to_nearest_atom) + '\n') 
                
               POSITION_FORCE_array[s*number_of_particles + o][6] = -1 # unbond so that we don't double record

     progess.updateAmount(s)
     print progess, "\r",

outputFile.close()

print
