
# coding: utf-8

# In[4]:


get_ipython().magic(u'pylab inline')
import mdtraj as md


# In[5]:

# Load parameters
top= md.load_prmtop("SYSTEM.top")


# In[12]:

# Select residues of ATP phosphate and Tyr
atp_phos = top.select("resSeq 288 and name PG") 
tyr = top.select("resSeq 53 and name OH")


# In[13]:

# Check correct atoms selected
atom = top.atom(int(atp_phos))   
print('''This is the %sth atom, atom type is %s.
Part of a %s residue.''' % ( atom.index, atom.name, atom.residue.name))
atom = top.atom(int(tyr))   
print('''This is the %sth atom, atom type is %s.
Part of a %s residue.''' % ( atom.index, atom.name, atom.residue.name))


# In[15]:

# Make an array containing each atom pair involved in the interaction
# For md.compute_distances need a (Any, 2) shape array.
phos_tyr = np.array([atp_phos,tyr]) 
phos_tyr_new = phos_tyr.reshape(1,2) 


# In[17]:

distance_list_1 = []
for chunk in md.iterload('TRAJ.dcd', chunk=100, top=top):
    distance_list_1.append(md.compute_distances(chunk, phos_tyr_new, periodic=True, opt=True))


# In[22]:

# Convert lists to arrays
distance_array_1 = np.array(distance_list_1)
distance_array_new_1 = np.reshape(distance_array_1, (10000,1))

# Convert to angstrom
phos_tyr_angstrom = distance_array_new_1 * 10


# In[23]:

# Ensure appropriate histogram range selected
print "RANGE OF VALUES"
print "min", np.amin(phos_tyr_angstrom), "max", np.amax(phos_tyr_angstrom)


# In[24]:

min_bin = 7.00
max_bin = 25.00

print "Using bin range", min_bin, "to", max_bin, "Type NO if this is not suitable based on max/min values above. Any other key to continue."
if raw_input() =="NO":
	print "Stopped. Edit range of histogram and run again."
	sys.exit()


# In[30]:

# Plot histogram and save as txt file
print "Histogram of ATP gamma phosphate to Tyrosine126 O(H) distance."
(n, bins) = np.histogram(phos_tyr_angstrom, bins = 100, range = (min_bin, max_bin), normed=True)
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
np.savetxt("TYR_OUTPUT/output.dat", data, fmt=['%d', '%.20f'])


# In[31]:

x = np.mean(distance_array_new_1)
print "The average distance between ATP gamma-phosphate and Tyr126 OH is %s nm." % (x)
x_ang = x * 10
x_ang_c = "%.3f" % x_ang
print "%s Angstroms." % (x_ang_c)


# In[ ]:



