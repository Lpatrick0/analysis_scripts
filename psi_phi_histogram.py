# Run from a directory containing folders organised as follows:
# /OUTPUT/PSI
# /OUTPUT/PHI

# Run as $ python psi_phi_histograms.py parameterfile.parm7 trajectoryfile.dcd

import scipy as sp
from scipy import stats
import numpy as np
import matplotlib.pyplot as plt
import mdtraj as md
import sys
from sys import argv
import math
#get_ipython().magic(u'pylab inline')

script, parmfile, trajfile = argv

# Load parameter/topology file
top = md.load_prmtop(parmfile)

# Iterates over trajectory loading in chunks to save memory
for chunk in md.iterload(trajfile, chunk=100, top=top):
    psi_list = md.compute_psi(chunk)
    phi_list = md.compute_phi(chunk)

# Psi angles histogram
for i in xrange(psi_list[1].shape[1]):
    dihedral_traj = psi_list[1][:,i]
    dihedral_traj_deg = (dihedral_traj * 180) / math.pi
    index = np.linspace(1,len(dihedral_traj), len(dihedral_traj)).astype('int')
    # Histogram
    (n, bins) = np.histogram(dihedral_traj_deg, bins = 300, range = (-180.00, 180.00), normed=True)
    bincentre = 0.5*(bins[1:]+bins[:-1])
    index = np.linspace(1, len(bincentre), num = len(bincentre), dtype = int)
    # Total amount to be split over empty bins only
    total_bin_addition = 0.001
    all_bins = len(bincentre)
    non_zero = np.count_nonzero(n)
    zero_bins = all_bins - non_zero
    bin_addition = total_bin_addition/float(zero_bins)
    # Adds the bin_addition amount into all zero-count bins
    for j in xrange(len(n)):
        if n[j] == 0.0:
            n[j] = bin_addition
    #normalise
    n = n / (sum(n))
    data = np.vstack((index, n)).T
    np.savetxt('OUTPUT/PSI/psi_hist_%d.dat' % (i+1), data, fmt=['%d', '%.30f'])

# Phi angles histogram
for i in xrange(phi_list[1].shape[1]):
    dihedral_traj = phi_list[1][:,i]
    dihedral_traj_deg = (dihedral_traj * 180) / math.pi
    index = np.linspace(1,len(dihedral_traj), len(dihedral_traj)).astype('int')
    (n, bins) = np.histogram(dihedral_traj_deg, bins = 300, range = (-180.00, 180.00), normed=True)
    bincentre = 0.5*(bins[1:]+bins[:-1])
    index = np.linspace(1, len(bincentre), num = len(bincentre), dtype = int)
    total_bin_addition = 0.001
    all_bins = len(bincentre)
    non_zero = np.count_nonzero(n)
    zero_bins = all_bins - non_zero
    bin_addition = total_bin_addition/float(zero_bins)
    for j in xrange(len(n)):
        if n[j] == 0.0:
            n[j] = bin_addition
    n = n / (sum(n))
    data = np.vstack((index, n)).T
    np.savetxt('OUTPUT/PHI/psi_hist_%d.dat' % (i+1), data, fmt=['%d', '%.30f'])
