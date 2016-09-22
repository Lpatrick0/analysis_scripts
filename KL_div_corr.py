# use pylab inline with ipython, but switch to import when converting to a script
# get_ipython().magic(u'pylab inline')
# -*- coding: utf-8 -*-
import scipy as sp
from scipy import stats
import numpy as np
import matplotlib.pyplot as plt

FILENAME_1 = raw_input("Reference ensemble filename-> ")
FILENAME_2 = raw_input("Target ensemble filename-> ")
FILENAME_3 = raw_input("Output file name->")

# Load two text files as np arrays, produced by 'distance_between_atoms_and_print_hist_data.py'.
P = np.loadtxt(FILENAME_1)
Q = np.loadtxt(FILENAME_2)

KL = sp.stats.entropy(pk=P, qk=Q, base=None)

print KL
np.savetxt(FILENAME_3, KL)
