#Determination of distance between a glutamic acid residue, located on the alpha-C helix,
#  and a key active site lysine residue which interacts with a phosphate of ATP.
#Atoms chosen to measure distance NZ of Lys and CD of Glu.

# coding: utf-8
# get_ipython().magic(u'pylab inline')
import scipy as sp
from scipy import stats
import numpy as np
import matplotlib.pyplot as plt
import mdtraj as md
import sys
from sys import argv
# use pylab inline with ipython, but switch to import when converting to a script

script, parmfile, trajfile, output = argv

# Load topology file
top= md.load_prmtop(parmfile)

#If required, can edit to select residues (resSeq) & atoms (using atom name).
res1 = top.select("resSeq 38 and name NZ")
res2 = top.select("resSeq 57 and name CD")

#To check using correct atoms
atom = top.atom(res1)
print('''First atom for distance measurement is the %sth atom, atom type is %s.
Part of a %s residue.''' % ( atom.index, atom.name, atom.residue.name))
atom = top.atom(res2)
print('''Second atom for distance measurement is the %sth atom, atom type is %s.
Part of a %s residue.''' % ( atom.index, atom.name, atom.residue.name))

# Make an array containing each atom pair involved in the interaction.
atom_pairs = np.array([res1,res2])    # make atom pairs

# For md.compute_distances need a (Any, 2) shape array.
new = atom_pairs.reshape(1,2)
print new

distance_list = []
#loads md trajectory in chunks (saves memory) with parameter file
for chunk in md.iterload(trajfile, chunk=100, top=top):
    distance_list.append(md.compute_distances(chunk, new, periodic=True, opt=True))

# Combine data from each chunk into 1 array
distance_array = np.array(distance_list)
distance_array_new = np.reshape(distance_array, (10000,1))

# Plot histogram
(n, bins) = np.histogram(distance_array_new, bins = 100, range = (0.00, 1.00), normed=True)

# Find centre of each bin (rather than edge of bin)
bincentre = 0.5*(bins[1:]+bins[:-1])

# Create an array the length of the number of bins
index = np.linspace(1, len(bincentre), num = len(bincentre), dtype = int)

# Add to empty bins only - amount added is equal to (total bin addition)/(number empty bins) - to allow calculation of KL for non-continuous distributions
total_bin_addition = 0.001

# Subtracts the number of non-empty bins from total bins to get the number of empty bins, toa llow the total_bin_addition to be divided
all_bins = len(bincentre)
non_zero = np.count_nonzero(n)
print "number of non-zero bins: ", non_zero
zero_bins = all_bins - non_zero
print "number of zero bins: ", zero_bins
bin_addition = total_bin_addition/float(zero_bins)
print "Value added to empty bins: ", bin_addition

# Adds the bin_addition amount into all zero-count bins
for i in xrange(len(n)):
    if n[i] == 0.0:
        n[i] = bin_addition

#normalise
n = n / (sum(n))

# Arrange index and frequency of each bin in array
data = np.vstack((index, n)).T
np.savetxt(output, data, fmt=['%d', '%.10f'])


x = np.mean(distance_array_new)
print "The average distance between NZ (Lys37) and CD (Glu56) is %s nm." % (x)
x_ang = x * 10
x_ang_c = "%.3f" % x_ang
print "%s Angstroms." % (x_ang_c)
