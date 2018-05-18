#!/Library/Frameworks/Python.framework/Versions/2.7/bin/python
# ----------------------------------------
# COMPUTE THE INVERSE NxN COVARIANCE MATRIX AND AVERAGE POSITIONS FROM A TRAJECTORY
# ----------------------------------------
# USAGE:

# ./compute_covar_avg.py config_file

# ----------------------------------------
# PREAMBLE:


import sys
import numpy as np
from numpy.linalg import *

zeros = np.zeros
sqrt = np.sqrt
sum = np.sum
eigen = np.linalg.eig
flush = sys.stdout.flush


# ----------------------------------------
# VARIABLE DECLARATION:


config_file = sys.argv[1]

necessary_parameters = ['wt_average_file','wt_covar_inv_file','mut_average_file','mut_res_num','res_min','res_max']
all_parameters = ['wt_average_file','wt_covar_inv_file','mut_average_file','mut_res_num','res_min','res_max']


# ----------------------------------------
# SUBROUTINES:

def ffprint(string):
    print '%s' %(string)
    flush()

def config_parser(config_file): # Function to take config file and create/fill the parameter dictionary
        for i in range(len(necessary_parameters)):
                parameters[necessary_parameters[i]] = ''

        # GRABBING PARAMETER VALUES FROM THE CONFIG FILE:
        execfile(config_file,parameters)
        for key, value in parameters.iteritems():
                if value == '':
                        print '%s has not been assigned a value. This variable is necessary for the script to run. Please declare this variable within the config file.' %(key)
                        sys.exit()



# ----------------------------------------
# MAIN:                                           
# ----------------------------------------
# CREATING PARAMETER DICTIONARY
parameters = {}
config_parser(config_file)


# get inverse covariance matrix (from wt simulation)
covar_inv = np.loadtxt(parameters['wt_covar_inv_file'])
# get wt average positions
wt_avg_pos = np.loadtxt(parameters['wt_average_file'])
# get mutant average positions
mut_avg_pos = np.loadtxt(parameters['mut_average_file'])


n_res = wt_avg_pos.shape[0]
lnKd = 0.0
for j in range(int(parameters['res_min'])-1,int(parameters['res_max'])):
	delta_r_j = mut_avg_pos[j,:] - wt_avg_pos[j,:]
	for i in range(n_res):
		delta_r_i = mut_avg_pos[i,:] - wt_avg_pos[i,:]
		lnKd += covar_inv[i,j] * np.dot(delta_r_j,delta_r_i)
# average Kd over observation points in the j loop
lnKd /= float( int(parameters['res_max']) - int(parameters['res_min']) + 1 )
# print Kd
print "ln(K_d):", lnKd






