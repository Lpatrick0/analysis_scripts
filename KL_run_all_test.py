# use pylab inline with ipython, but switch to import when converting to a script
# get_ipython().magic(u'pylab inline')
# -*- coding: utf-8 -*-
import scipy as sp
from scipy import stats
import numpy as np
import matplotlib.pyplot as plt
from sys import argv

script, reference, target, output = argv

# Load two text files as np arrays
P = np.loadtxt(reference)
Q = np.loadtxt(target)

KL = sp.stats.entropy(pk=P, qk=Q, base=None)

print KL
np.savetxt(output, KL)
