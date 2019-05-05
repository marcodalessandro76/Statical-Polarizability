# Before calling this script (for istance in the launcher slurm script), perform the operatotions:
# export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
# export BIGDFT_MPIRUN='mpirun -np $value$'


from BigDFT import Datasets as D, Calculators as C, Inputfiles as I, Logfiles as lf
from BigDFT.Database import Molecules
from futile.Utils import write
import numpy as np
import matplotlib.pyplot as plt
import os, sys, yaml
sys.path.insert(0,'../')
import StatPol as SP, Routines as R

basepath = os.getcwd()

nsp_dataset = yaml.load(open('nsp_dataset.yaml')) #load the dataset with the list of gs_study
#nsp_dataset = yaml.load(open('nsp_results.yaml')) #load the dataset with the results already computed

# For testing and computation purposes we split the dataset into several subset made of
# 10 molecules each (molecules are sorted in alphabetical order before the splitting). So we have
# 8 subsets (the last one with 5 molecules)

molecules = nsp_dataset.keys()
molecules.sort()
subset_length = 10
subset = {}
for ind in range(int(len(molecules)/subset_length)+1):
    subset[ind] = []
for ind,mol in enumerate(molecules):
    subset[int(ind/subset_length)].append(mol)

reduced_study_set = [('lda_pt','hgh_k')]#,('pbe0','hgh_k')]
#('lda_pw','hgh_k'),('pbe','hgh_k'),('pbe','nlcc_aw'),('pbe','nlcc_ss')]

reload(R)
# mpi and omp are set from above exporting the asscoiated variables
code=C.SystemCalculator(skip=True,verbose=False)

for mol in ['CO']# subset[0]:
    data = nsp_dataset[mol]
    data['results'] = {}
    for study in data['study']:
        if study in reduced_study_set:
            data['results'][study] = {}
            data['results'][study] = R.nsp_workflow(term_verb=True,molecule=mol,study=study,code=code,data_folder='data_hpc')
            print ''
