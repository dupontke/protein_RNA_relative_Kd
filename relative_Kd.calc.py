#!/Library/Frameworks/Python.framework/Versions/2.7/bin/python
# ----------------------------------------
# USAGE:

# ./relative_Kd.calc.py config_file

# ----------------------------------------
# PREAMBLE:


import sys
import numpy as np
from numpy.linalg import *
import MDAnalysis
import MDAnalysis.analysis.align

zeros = np.zeros
sqrt = np.sqrt
sum = np.sum
eigen = np.linalg.eig
flush = sys.stdout.flush


# ----------------------------------------
# VARIABLE DECLARATION:


config_file = sys.argv[1]

necessary_parameters = ['pdb_file','traj_loc','start','end','average_pdb']
all_parameters = ['pdb_file','traj_loc','start','end','average_pdb','alignment','covar_selection','coarseness','fine_grain_selection','cartesian_correlation_filename','cartesian_average_filename','cartesian_variance_filename','cartesian_covariance_filename','distance_correlation_filename','distance_variance_filename','distance_covariance_filename','functionalize_distance_correlation_bool','functionalized_distance_correlation_filename','PCA_bool','PCA_eigenvalues_filename','PCA_eigenvectors_filename','summary_bool','summary_filename']


# ----------------------------------------
# SUBROUTINES:

def ffprint(string):
    print '%s' %(string)
    flush()

def config_parser(config_file): # Function to take config file and create/fill the parameter dictionary
        for i in range(len(necessary_parameters)):
                parameters[necessary_parameters[i]] = ''

        # SETTING DEFAULT PARAMETERS FOR OPTIONAL PARAMETERS:
        parameters['alignment'] = 'protein'
        parameters['covar_selection'] = 'protein'
        parameters['coarseness'] = 'COM'
        parameters['fine_grain_selection'] = None
        parameters['cartesian_correlation_filename'] = 'cartesian_correlation.dat'
        parameters['cartesian_average_filename'] = 'cartesian_average.dat'
        parameters['cartesian_variance_filename'] = 'cartesian_variance.dat'
        parameters['cartesian_covariance_filename'] = 'cartesian_covariance.dat'
        parameters['distance_correlation_filename'] = 'distance_correlation.dat'
        parameters['distance_variance_filename'] = 'distance_variance.dat'
        parameters['distance_covariance_filename'] = 'distance_covariance.dat'
        parameters['functionalize_distance_correlation_bool'] = False
        parameters['functionalized_distance_correlation_filename'] = 'functionalized_dist_covar.dat'
        parameters['PCA_bool'] = False
        parameters['PCA_eigenvalues_filename'] = 'PCA_eigenvalues_cartesian_covariance.dat'
        parameters['PCA_eigenvectors_filename'] = 'PCA_eigenvectors_cartesian_covariance.dat'
        parameters['summary_bool'] = True
        parameters['summary_filename'] = 'water_retention_analysis.summary'

        # GRABBING PARAMETER VALUES FROM THE CONFIG FILE:
        execfile(config_file,parameters)
        for key, value in parameters.iteritems():
                if value == '':
                        print '%s has not been assigned a value. This variable is necessary for the script to run. Please declare this variable within the config file.' %(key)
                        sys.exit()

        if parameters['coarseness'] not in ['COM','Atomic']:
                print "coarseness parameter does not match an acceptable value. Viable values are 'COM' and 'Atomic'. Killing job."
                sys.exit()

def summary(filename):
        with open(filename,'w') as W:
                W.write('Using MDAnalysis version: %s\n' %(MDAnalysis.version.__version__))
                W.write('To recreate this analysis, run this line:\n')
                for i in range(len(sys.argv)):
                        W.write('%s ' %(sys.argv[i]))
                W.write('\n\nParameters used:\n')
                for i in all_parameters:
                        W.write('%s = %s \n' %(i,parameters[i]))
                W.write('\n\n')


# ----------------------------------------
# MAIN:                                           
# ----------------------------------------
# CREATING PARAMETER DICTIONARY
parameters = {}
config_parser(config_file)

start = int(parameters['start'])
end = int(parameters['end'])

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











#------------------------------------------
# OLLD script... need to convert
#------------------------------------------

pdb = sys.argv[1]                   # point to a pdb or prmtop or psf file (untested for both prmtop and psf files)
traj_loc = sys.argv[2]              # point to the location of the trajectory files
start = int(sys.argv[3])            # integer describing
end = int(sys.argv[4])
system = sys.argv[5]


important = 'protein and resname GTP and resname SAH and resname MG'
nSel = len(sel)

# ----------------------------------------
# FUNCTIONS:

def ffprint(string):
    print '%s' %(string)
    flush()

# ----------------------------------------
# MAIN PROGRAM:

u = MDAnalysis.Universe(pdb)
u_important = u.select_atoms(important)

nRes = len(u_important.residues)
ffprint(nRes)

# make an atom selection to compute covariance matrix distance
u_sel = []
for i in range(nSel):
    selection0 = sel[i][0]
    selection1 = sel[i][1]
    u_sel.append([u.select_atoms(selection0),u.select_atoms(selection1)])
    ffprint('%s atoms found in %s' %(u_sel[i][0].n_atoms, u_sel[i][0].residues))
    ffprint('%s atoms found in %s' %(u_sel[i][1].n_atoms, u_sel[i][1].residues))





out1 = open('%03d.%03d.com_distance.dat' %(int(sys.argv[3]),end),'w')

count = 0
nSteps = 0
while start <= end:
    ffprint('Loading trajectory %s' %(start))
    u.load_new('%sproduction.%s/production.%s.dcd' %(traj_loc,start,start))
    nSteps += len(u.trajectory)

    for ts in u.trajectory:
        if ts.frame%1000 == 0:
            ffprint('Working on timestep %d of trajectory %d' %(ts.frame, start))

        for i in range(nSel):

            out1.write('%10.6f    ' %(dist))
        out1.write('\n')
    start +=1

out1.close()
ffprint(nSteps)
