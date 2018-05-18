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
import MDAnalysis
from MDAnalysis.analysis import align

zeros = np.zeros
sqrt = np.sqrt
sum = np.sum
eigen = np.linalg.eig
flush = sys.stdout.flush


# ----------------------------------------
# VARIABLE DECLARATION:


config_file = sys.argv[1]

necessary_parameters = ['top_file','traj_file', 'ref_pdb_file','start','end','average_out','covar_out','inv_covar_out']
all_parameters = ['top_file','traj_file', 'ref_pdb_file','start','end','average_out','covar_out','inv_covar_out','alignment','covar_selection']


# ----------------------------------------
# SUBROUTINES:

def ffprint(string):
    print '%s' %(string)
    flush()

def config_parser(config_file): # Function to take config file and create/fill the parameter dictionary
        for i in range(len(necessary_parameters)):
                parameters[necessary_parameters[i]] = ''

        # SETTING DEFAULT PARAMETERS FOR OPTIONAL PARAMETERS:
        parameters['alignment'] = 'name CA'
        parameters['covar_selection'] = 'name CA'

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

#start = int(parameters['start'])
#end = int(parameters['end'])

# ----------------------------------------
# OUTLINE FOR CALCULATIONS:
# 1. reference point for force calculations is the mutation site
# 2. electrostatic force is F(i) is generated on atom(i) by the charges of all the atoms on the sidechain of the mutated residue.
#        a. the force varies before and after the mutation dF(i) = F(i)wt - F(i)mut
#     I. so in other words I need to calculate the electrostatic force for all atoms with the reference point at the mutated residue. I need to do this for the wt and the mutation in question. Then take the difference between the force for that atom atom pair -> aka covariance  

# what ensemble are they averaging over for the covariance matrix??
#---> the code up covariance via a dot product to get NxN instead of 3Nx3N. then calculate covariance. need to know for eq 5 is the inverted covariance and the arrary of average ca position for WT and mutant.
## NxN covariance and inverted covariance and the average positions. all we need for the mutants are the average positions.
# read in the inverted covariance of wt, average positiong for wt and mutants  





# Initiate MD Analysis universe

u = MDAnalysis.Universe(parameters['top_file'], parameters['traj_file'])
u_align = u.select_atoms(parameters['alignment'])

ref = MDAnalysis.Universe(parameters['ref_pdb_file'])
ref_align = ref.select_atoms(parameters['alignment'])

avg_pos = np.zeros((u_align.n_atoms,3),dtype=float)
covar = np.zeros((u_align.n_atoms,u_align.n_atoms),dtype=float)

print u.trajectory.n_frames

for ts in u.trajectory:
	# align frame to reference
	align.alignto(u,ref,select=parameters['alignment'])
	# add to average
	avg_pos += u_align.positions
	# add to covar
	covar += np.dot(u_align.positions,u_align.positions.T)

# finalize average
avg_pos /= float(u.trajectory.n_frames)

# finalize covariance
covar = covar/float(u.trajectory.n_frames) - np.dot(avg_pos,avg_pos.T)

# print out averages
np.savetxt(parameters['average_out'], avg_pos)

# print out covariance
np.savetxt(parameters['covar_out'],covar)

# 
