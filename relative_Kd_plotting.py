#!/Library/Frameworks/Python.framework/Versions/2.7/bin/python
# USAGE:

# PREAMBLE:

import numpy as np
import sys
import os
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib as mpl
from plotting_functions import *
from sel_list import *

zeros = np.zeros

dat = sys.argv[1]  
system = sys.argv[2]

nSel = len(sel)

bin_size = 0.01
k = 0.001987 # Kcal K^-1 mol^-1
T = 300. # K
kT = k*T
boltz = 2*kT
four_pi = 4*np.pi


# ----------------------------------------
# MAIN PROGRAM:
# ----------------------------------------

# Load in data_file into a numpy array
datalist = np.loadtxt(dat)

nSteps = len(datalist[:,0])          
print 'Number of selections: %d, number of steps: %d' %(nSel,nSteps)

time = np.zeros(nSteps)
for i in range(nSteps):
	time[i] = i*0.002		# units of time in ns; each frame is separated by 0.002 ns 

