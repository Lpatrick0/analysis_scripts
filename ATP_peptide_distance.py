# Run from a folder with directory structure *current folder*/ATP_OUTPUT/


# coding: utf-8
# get_ipython().magic(u'pylab inline')
import scipy as sp
from scipy import stats
import numpy as np
import matplotlib.pyplot as plt
import mdtraj as md
import sys
from sys import argv
#use pylab inline with ipython, but switch to import when converting to a script
import mdtraj as md

# Load parameters
top= md.load_prmtop("complex_solv.parm7")

# Load trajectory
traj = md.load('traj000000002.dcd', top=top)
print traj

# Select residues of ATP phosphate and Peptide Thr
atp_phos = top.select("resSeq 288 and name PG") 
pep_thr = top.select("resSeq 293 and name OG1")

# Check correct atoms selected
atom = top.atom(int(atp_phos))   
print('''This is the %sth atom, atom type is %s.
Part of a %s residue.''' % ( atom.index, atom.name, atom.residue.name))
atom = top.atom(int(pep_thr))   
print('''This is the %sth atom, atom type is %s.
Part of a %s residue.''' % ( atom.index, atom.name, atom.residue.name))

# Make an array containing each atom pair involved in the interaction
# For md.compute_distances need a (Any, 2) shape array.
phos_thr = np.array([atp_phos,pep_thr]) 
phos_thr_new = phos_thr.reshape(1,2) 

# Compute distances
phos_thr_distance = md.compute_distances(traj, phos_thr_new, periodic=True, opt=True)
# Convert to angstrom
phos_thr_angstrom = phos_thr_distance * 10

# Ensure appropriate histogram range selected
print "RANGE OF VALUES"
print "min", np.amin(phos_thr_angstrom), "max", np.amax(phos_thr_angstrom)

min_bin = 1.00
max_bin = 15.00

print "Using bin range", min_bin, "to", max_bin, "Type NO if this is not suitable based on max/min values above. Any other key to continue."
if raw_input() =="NO":
	print "Stopped. Edit range of histogram and run again."
	sys.exit()

# Plot histogram and save as txt file
print "Histogram of ATP gamma phosphate to Peptide Threonine distance."
(n, bins) = np.histogram(phos_thr_angstrom, bins = 100, range = (min_bin, max_bin), normed=True)
n = n / (sum(n))
bincentre = 0.5*(bins[1:]+bins[:-1])
index = np.linspace(1, len(bincentre), num = len(bincentre), dtype = int)
# Add to empty bins only - amount added is equal to 1/(number empty bins) - to allow calculation of KL
total_bin_addition = 0.000001
all_bins = len(bincentre)
# To count the number of populated and non populated bins, to allow dividion of the total bin addition
non_zero = np.count_nonzero(n)
print "Number of populated bins:", non_zero
zero_bins = all_bins - non_zero
print "Number of zero bins:", zero_bins
bin_addition = total_bin_addition/float(zero_bins)
print "Amount added to empty bins: ", bin_addition
for i in xrange(len(n)):
    if n[i]==0.0:
        n[i] = bin_addition
data = np.vstack((index, n)).T
np.savetxt("ATP_OUTPUT/output.dat", data, fmt=['%d', '%.20f'])

x = np.mean(phos_thr_angstrom)
print "The average distance between PG (ATP) and OG1 (THR-PEPTIDE) is %s Angstroms." % (x)
