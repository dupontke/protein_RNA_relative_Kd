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

necessary_parameters = ['pdb_file','ref_pdb_file','traj_loc','start','end','alignment','covar_selection','average_out','covar_out','inv_covar_out']
all_parameters = ['pdb_file','ref_pdb_file','traj_loc','start','end','alignment','covar_selection','average_out','covar_out','inv_covar_out','alignment','covar_selection','write_summary','summary_filename']


# ----------------------------------------
# SUBROUTINES:

def ffprint(string):
    print '%s' %(string)
    flush()

def config_parser(config_file): # Function to take config file and create/fill the parameter dictionary
        for i in range(len(necessary_parameters)):
                parameters[necessary_parameters[i]] = ''

        # SETTING DEFAULT PARAMETERS FOR OPTIONAL PARAMETERS:
        parameters['write_summary'] = False
        parameters['summary_filename'] = 'compute_covar_avg.summary'

        # GRABBING PARAMETER VALUES FROM THE CONFIG FILE:
        execfile(config_file,parameters)
        for key, value in parameters.iteritems():
                if value == '':
                        print '%s has not been assigned a value. This variable is necessary for the script to run. Please declare this variable within the config file.' %(key)
                        sys.exit()

def summary():
    with open('%s' %(parameters['summary_filename']),'w') as f:
        f.write('Using MDAnalysis version: %s\n' %(MDAnalysis.version.__version__))
        f.write('To recreate this analysis, run this line:\n')
        for i in range(len(sys.argv)):
            f.write('%s' %(sys.argv[i]))
        f.write('\n\n')
        f.write('Parameters used:\n')
        for i in all_parameters:
            f.write('%s = %s \n' %(i,parameters[i]))
        f.write('\n\n')



# ----------------------------------------
# MAIN:                                           
# ----------------------------------------
# CREATING PARAMETER DICTIONARY
parameters = {}
config_parser(config_file)

# ----------------------------------------
# Initiate MD Analysis universe
ref = MDAnalysis.Universe(parameters['ref_pdb_file'])
ref_align = ref.select_atoms(parameters['alignment'])

u = MDAnalysis.Universe(parameters['pdb_file'])
u_align = u.select_atoms(parameters['alignment'])

# ----------------------------------------
# ARRAY DECLARATION
avg_pos = np.zeros((u_align.n_atoms,3),dtype=float)
covar = np.zeros((u_align.n_atoms,u_align.n_atoms),dtype=float)

# ----------------------------------------
# TRAJECTORY ANALYSIS
nSteps = 0
start = int(parameters['start'])
end = int(parameters['end'])
while start <= end:
    ffprint('Loading trajectory %s' %(start))
    u.load_new('%sproduction.%s/production.%s.dcd' %(parameters['traj_loc'],start,start))
    nSteps += len(u.trajectory)
    print nSteps
    # Loop through trajectory
    for ts in u.trajectory:
        # align frame to reference
        align.alignto(u,ref,select=parameters['alignment'])
        # add to average
        avg_pos += u_align.positions
        # add to covar
        covar += np.dot(u_align.positions,u_align.positions.T)
    start += 1

# finalize average
avg_pos /= nSteps
#avg_pos /= float(u.trajectory.n_frames)

# finalize covariance
covar = covar/nSteps - np.dot(avg_pos,avg_pos.T)
#covar = covar/float(u.trajectory.n_frames) - np.dot(avg_pos,avg_pos.T)

# get inverse of covariance matrix
inv_covar = np.linalg.inv(covar)



# print out averages
np.savetxt(parameters['average_out'], avg_pos)

# print out covariance
np.savetxt(parameters['covar_out'],covar)

# print out inverse of covariance matrix
np.savetxt(parameters['inv_covar_out'],inv_covar)


if parameters['write_summary']:
    summary()
