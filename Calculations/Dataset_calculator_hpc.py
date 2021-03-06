# Before calling this script (for istance in the launcher slurm script), perform the operations:
# export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
# export BIGDFT_MPIRUN='mpirun -np $value$'

# Run the script by adding one (or more) integers between 0 and 7 that specify the number of the
# subsets of molecules that are computed

from BigDFT import Datasets as D, Calculators as C, Inputfiles as I, Logfiles as lf
from BigDFT.Database import Molecules
from futile.Utils import write
import numpy as np
import os, sys, yaml
sys.path.insert(0,'../')
import StatPol as SP, Routines as R

basepath = os.getcwd()

nsp_dataset = yaml.load(open('nsp_dataset.yaml')) #load the dataset with the list of gs_study
#nsp_dataset = yaml.load(open('nsp_results.yaml')) #load the dataset with the results already computed

# For testing and computation purposes we split the dataset into several subset made of
# 10 molecules each (molecules are sorted in alphabetical order before the splitting). So we have
# 8 subsets (the last one with 4 molecules)

molecules = nsp_dataset.keys()
molecules.sort()
subset_length = 10
subset = {}
for ind in range(int(len(molecules)/subset_length)+1):
    subset[ind] = []
for ind,mol in enumerate(molecules):
    subset[int(ind/subset_length)].append(mol)

calc = []
for arg in sys.argv[1:]:
    calc += subset[int(arg)]
print 'Compute molecules'
print calc

reduced_study_set = [('lda_pt','hgh_k'),('lda_pw','hgh_k'),('pbe','hgh_k'),('pbe','nlcc_aw'),\
                     ('pbe','nlcc_ss'),('pbe0','hgh_k')]

# mpi and omp are set from above exporting the asscoiated variables
code=C.SystemCalculator(skip=True,verbose=True)

for mol in calc:
    data = nsp_dataset[mol]
    data['results'] = {}
    for study in data['study']:
        if study in reduced_study_set:
            data['results'][study] = {}
            data['results'][study] = R.nsp_workflow(alpha_conv=1.0e-3,term_verb=True,molecule=mol,study=study,code=code,data_folder='Data')
            print ''

# Save the dataset as yaml file
#with open('nsp_results_tol=1em3.yaml', 'w') as outfile:
#    yaml.dump(nsp_dataset, outfile, default_flow_style=False)
