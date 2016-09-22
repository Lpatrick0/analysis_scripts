# -*- coding: utf-8 -*-

import numpy as np

KL_list = []
for i in xrange(0, 285):
    KL = np.loadtxt("KL_all_angles/output/psi_%s.dat" % i)
    KL_list.append(KL[1])

KL_list_2dp = [ '%.2f' % elem for elem in KL_list ]
print KL_list_2dp

np.savetxt('1_list_KL_psi/KL_psi.dat', KL_list_2dp, fmt="%s")
