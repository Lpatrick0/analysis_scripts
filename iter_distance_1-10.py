# coding: utf-8

#Determination of distance between a glutamic acid residue, located on the alpha-C helix,
#  and a key active site lysine residue which interacts with a phosphate of ATP.
#Atoms chosen to measure distance NZ of Lys and CD of Glu.

# Should run from a folder containing 10 trajectories (named as traj000000001.dcd, traj000000002.dcd, etc) and parameter file. Run as $ python *this_script*.py paramfile.parm7 output_filename.dat


# get_ipython().magic(u'pylab inline')
import scipy as sp
from scipy import stats
import numpy as np
import matplotlib.pyplot as plt
import mdtraj as md
import sys
from sys import argv
#use pylab inline with ipython, but switch to import when converting to a script

script, parmfile, output = argv

# Load topology / parameter file
top= md.load_prmtop(parmfile)

residue1_input = raw_input("First residue number: ")
residue1_input_name = raw_input("First residue name: ")
residue2_input = raw_input("Second residue number: ")
residue2_input_name = raw_input("Second residue name: ")

first = "resSeq " + residue1_input + " and name " + residue1_input_name
second = "resSeq " + residue2_input + " and name " + residue2_input_name

res1 = top.select(first)
res2 = top.select(second)

# Check correct atoms selected
atom = top.atom(res1)
print('''This is the %sth atom, atom type is %s.
Part of a %s residue.''' % (atom.index, atom.name, atom.residue.name))
atom = top.atom(res2)
print('''This is the %sth atom, atom type is %s.
Part of a %s residue.''' % (atom.index, atom.name, atom.residue.name))

# Make an array containing each atom pair involved in the interaction.
atom_pairs = np.array([res1,res2])
# For md.compute_distances need a (Any, 2) shape array.
new = atom_pairs.reshape(1,2)
print new

# Loads trajectory in chunks and compute.distances adds to each list, for each trajectory (ten 100ns trajectories, each started from restart file from the last)
distance_list_1 = []
for chunk in md.iterload('traj000000001.dcd', chunk=100, top=top):
    distance_list_1.append(md.compute_distances(chunk, new, periodic=True, opt=True))

distance_list_2 = []
for chunk in md.iterload('traj000000002.dcd', chunk=100, top=top):
    distance_list_2.append(md.compute_distances(chunk, new, periodic=True, opt=True))

distance_list_3 = []
for chunk in md.iterload('traj000000003.dcd', chunk=100, top=top):
    distance_list_3.append(md.compute_distances(chunk, new, periodic=True, opt=True))

distance_list_4 = []
for chunk in md.iterload('traj000000004.dcd', chunk=100, top=top):
    distance_list_4.append(md.compute_distances(chunk, new, periodic=True, opt=True))

distance_list_5 = []
for chunk in md.iterload('traj000000005.dcd', chunk=100, top=top):
    distance_list_5.append(md.compute_distances(chunk, new, periodic=True, opt=True))

distance_list_6 = []
for chunk in md.iterload('traj000000006.dcd', chunk=100, top=top):
    distance_list_6.append(md.compute_distances(chunk, new, periodic=True, opt=True))

distance_list_7 = []
for chunk in md.iterload('traj000000007.dcd', chunk=100, top=top):
    distance_list_7.append(md.compute_distances(chunk, new, periodic=True, opt=True))

distance_list_8 = []
for chunk in md.iterload('traj000000008.dcd', chunk=100, top=top):
    distance_list_8.append(md.compute_distances(chunk, new, periodic=True, opt=True))

distance_list_9 = []
for chunk in md.iterload('traj000000009.dcd', chunk=100, top=top):
    distance_list_9.append(md.compute_distances(chunk, new, periodic=True, opt=True))

distance_list_10 = []
for chunk in md.iterload('traj000000010.dcd', chunk=100, top=top):
    distance_list_10.append(md.compute_distances(chunk, new, periodic=True, opt=True))

# Each list will contain 100 "sub" lists, as each trajectory split into 100 chunks (each with 100 snapshots).
distance_array_1 = np.array(distance_list_1)
distance_array_new_1 = np.reshape(distance_array_1, (10000,1))

distance_array_2 = np.array(distance_list_2)
distance_array_new_2 = np.reshape(distance_array_2, (10000,1))

distance_array_3 = np.array(distance_list_3)
distance_array_new_3 = np.reshape(distance_array_3, (10000,1))

distance_array_4 = np.array(distance_list_4)
distance_array_new_4 = np.reshape(distance_array_4, (10000,1))

distance_array_5 = np.array(distance_list_5)
distance_array_new_5 = np.reshape(distance_array_5, (10000,1))

distance_array_6 = np.array(distance_list_6)
distance_array_new_6 = np.reshape(distance_array_6, (10000,1))

distance_array_7 = np.array(distance_list_7)
distance_array_new_7 = np.reshape(distance_array_7, (10000,1))

distance_array_8 = np.array(distance_list_8)
distance_array_new_8 = np.reshape(distance_array_8, (10000,1))

distance_array_9 = np.array(distance_list_9)
distance_array_new_9 = np.reshape(distance_array_9, (10000,1))

distance_array_10 = np.array(distance_list_10)
distance_array_new_10 = np.reshape(distance_array_10, (10000,1))

# Combine all the arrays into one big array
mega_array = numpy.concatenate((distance_array_new_1, distance_array_new_2, distance_array_new_3, distance_array_new_4, distance_array_new_5, distance_array_new_6, distance_array_new_7, distance_array_new_8, distance_array_new_9, distance_array_new_10), axis=0)

# Bin data into histogram
(n, bins) = np.histogram(mega_array, bins = 300, range = (0.00, 1.00), normed=True)

#normalise
n = n / (sum(n))

# Find centre of each bin.
bincentre = 0.5*(bins[1:]+bins[:-1])
#print bincentre

# Create an array the length of the number of bins
index = np.linspace(1, len(bincentre), num = len(bincentre), dtype = int)

# Add to empty bins only - amount added is equal to 1/(number empty bins) - to allow calculation of KL
total_bin_addition = 0.000001
all_bins = len(bincentre)
# To count the number of populated and non populated bins, to allow dividion of the total bin addition
non_zero = numpy.count_nonzero(n)
print "Number of populated bins:", non_zero
zero_bins = all_bins - non_zero
print "Number of zero bins:", zero_bins
bin_addition = total_bin_addition/float(zero_bins)
print "Amount added to empty bins: ", bin_addition
for i in xrange(len(n)):
    if n[i]==0.0:
        n[i] = bin_addition

# Arrange index and frequency of each bin in array
data = np.vstack((index, n)).T
print data
np.savetxt(output, data, fmt=['%d', '%.20f'])

x = np.mean(mega_array)
print "The average distance between NZ (Lys37) and CD (Glu56) is %s nm." % (x)
x_ang = x * 10
x_ang_c = "%.3f" % x_ang
print "%s Angstroms." % (x_ang_c)
