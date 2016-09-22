# coding: utf-8

import scipy as sp
from scipy import stats
import numpy as np
import matplotlib.pyplot as plt
import mdtraj as md
import sys
from sys import argv
import math

script, parmfile, trajfile = argv

top = md.load_prmtop(parmfile)
traj = md.load(trajfile, top=top)
print(traj)

psi_list = md.compute_psi(traj)
#phi_list = md.compute_phi(traj)
#chi1_list = md.compute_chi1(traj)
#chi2_list = md.compute_chi2(traj)

for i in xrange(psi_list[1].shape[1]):
    dihedral_traj = psi_list[1][:,i]
    dihedral_traj_deg = (dihedral_traj * 180) / math.pi
    index = np.linspace(1,len(dihedral_traj), len(dihedral_traj)).astype('int')
    (n, bins, patches) = plt.hist(dihedral_traj_deg, bins = 300, range = (-180.00, 180.00), normed=True)
    # data = np.vstack((index,dihedral_traj_deg)).T
    #    np.savetxt('dihedrals/psi_%d.xvg' %i, data, fmt=['%d', '%.4f'])
    # Add a uniform prior to whole distribution to allow calculation of KL
    n += ((np.mean(n))/100)
    #print n
    n = n / (sum(n))
    #print n
    # Find centre of each bin.
    bincentre = 0.5*(bins[1:]+bins[:-1])
    #print bincentre
    # Create an array the length of the number of bins
    index = np.linspace(1, len(bincentre), num = len(bincentre), dtype = int)
    #print ind
    # Arrange index and frequency of each bin in array
    data = np.vstack((index, n)).T
    np.savetxt('files/psi_hist_%d.dat' % i, data, fmt=['%d', '%.10f'])

#for i in xrange(phi_list[1].shape[1]):
#    dihedral_traj = phi_list[1][:,i]
#    dihedral_traj_deg = (dihedral_traj * 180) / math.pi
#    index = np.linspace(1,len(dihedral_traj), len(dihedral_traj)).astype('int')
#    (n, bins, patches) = plt.hist(dihedral_traj_deg, bins = 300, range = (-180.00, 180.00), normed=True)
#    n += ((np.mean(n))/100)
#    n = n / (sum(n))
#    bincentre = 0.5*(bins[1:]+bins[:-1])
#    index = np.linspace(1, len(bincentre), num = len(bincentre), dtype = int)
#    data = np.vstack((index, n)).T
#    np.savetxt('new/phi_hist_%d.dat' % i, data, fmt=['%d', '%.10f'])

#for i in xrange(chi1_list[1].shape[1]):
#    dihedral_traj = chi1_list[1][:,i]
#    dihedral_traj_deg = (dihedral_traj * 180) / math.pi
#    index = np.linspace(1,len(dihedral_traj), len(dihedral_traj)).astype('int')
#    (n, bins, patches) = plt.hist(dihedral_traj_deg, bins = 300, range = (-180.00, 180.00), normed=True)
#    n += ((np.mean(n))/100)
#    n = n / (sum(n))
#    bincentre = 0.5*(bins[1:]+bins[:-1])
#    index = np.linspace(1, len(bincentre), num = len(bincentre), dtype = int)
#    data = np.vstack((index, n)).T
#    np.savetxt('new/chi1_hist_%d.dat' % i, data, fmt=['%d', '%.10f'])

#for i in xrange(chi2_list[1].shape[1]):
#    dihedral_traj = chi2_list[1][:,i]
#    dihedral_traj_deg = (dihedral_traj * 180) / math.pi
#    index = np.linspace(1,len(dihedral_traj), len(dihedral_traj)).astype('int')
#    (n, bins, patches) = plt.hist(dihedral_traj_deg, bins = 300, range = (-180.00, 180.00), normed=True)
#    n += ((np.mean(n))/100)
#    n = n / (sum(n))
#    bincentre = 0.5*(bins[1:]+bins[:-1])
#    index = np.linspace(1, len(bincentre), num = len(bincentre), dtype = int)
#    data = np.vstack((index, n)).T
#    np.savetxt('new/chi2_hist_%d.dat' % i, data, fmt=['%d', '%.10f'])
