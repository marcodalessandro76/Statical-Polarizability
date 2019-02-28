from BigDFT import Datasets as D, Calculators as C, Inputfiles as I, Logfiles as lf
from BigDFT.Database import Molecules
from futile.Utils import write
import numpy as np
import os, sys
sys.path.insert(0,'../')
import StatPol as SP, workflow_hpc as w

molecule_set = ['CO','N2','H2O','Mg','Mg2','Ar','HCl','CH4','HF','NaCl','Ne','He']
xc_set = ['lda'] #['lda','pbe']
psp_set = ['hgh-k']

for molecule in molecule_set:
    for xc in xc_set:
        for psp in psp_set:
            w.single_study_workflow(term_verb=True,molecule=molecule,xc=xc,psp=psp)
