from BigDFT.Database import Molecules
from BigDFT import Datasets as D, Calculators as C, Inputfiles as I, Logfiles as lf
import StatPol as SP
import os, sys, numpy as np
from shutil import copyfile

AuToA = 0.5291772085**3

def get_atoms(mol):
    """
    Return a list of string with all the atoms of a molecule. For same molecules with
    a non standard name-scheme the attribution of the atoms is done by hand
    """
    if mol == 'SO-trip': return(['S','O'])
    if mol == 'CH2-t': return(['C','H'])
    if mol == 'FH-OH': return(['F','H','O'])
    if mol == 'H2O-Li': return(['H','O','Li'])

    atoms = []
    # remove numbers from str
    for c in mol:
        if c.isdigit():
            mol = mol.replace(c,'')
    # split elements to identify the atoms
    while True:
        if len(mol)==1:
            atoms.append(mol)
            break
        if len(mol)==2 and mol[1].islower():
            atoms.append(mol)
            break
        if mol[1].isupper():
            atoms.append(mol[:1])
            mol = mol[1:]
            continue
        if mol[1].islower():
            atoms.append(mol[:2])
            mol = mol[2:]
            continue
    # remove duplicates (and lose ordering)
    atoms = list(set(atoms))
    return atoms

def xc_lda_pt(inp,molecule,datadir):
    inp.set_xc(1)

def xc_lda_pw(inp,molecule,datadir):
    inp.set_xc('-001012')
    for atom in get_atoms(molecule):
        key = 'psppar.'+atom
        inp[key]={'Pseudopotential XC': 1}

def xc_pbe(inp,molecule,datadir):
    inp.set_xc(11)

def xc_pbe0(inp,molecule,datadir):
    inp.set_xc('PBE0')
    for atom in get_atoms(molecule):
        key = 'psppar.'+atom
        inp[key]={'Pseudopotential XC': 11}

basepath = os.getcwd()
psp_nlcc_aw_path = basepath+'/psp_nlcc_aw'
psp_nlcc_ss_path = basepath+'/psp_nlcc_ss'

#print basepath

def xc_pbe_nlcc_aw(inp,molecule,datadir):
    inp.set_xc(11)
    for atom in get_atoms(molecule):
        src = psp_nlcc_aw_path+'/psppar.'+atom
        dst = basepath+'/'+datadir+'/'+molecule+'/pbe-nlcc_aw/psppar.'+atom
        #dst = basepath+'/Data/'+molecule+'/pbe-nlcc_aw-test/psppar.'+atom
        copyfile(src,dst)

def xc_pbe_nlcc_ss(inp,molecule,datadir):
    inp.set_xc(11)
    for atom in get_atoms(molecule):
        src = psp_nlcc_ss_path+'/psppar.'+atom
        dst = basepath+'/'+datadir+'/'+molecule+'/pbe-nlcc_ss/psppar.'+atom
        #dst = basepath+'/Data/'+molecule+'/pbe-nlcc_ss-test/psppar.'+atom
        copyfile(src,dst)

set_xc = {('lda_pt','hgh_k') : xc_lda_pt,('lda_pw','hgh_k') : xc_lda_pw, \
         ('pbe','hgh_k') : xc_pbe, ('pbe0','hgh_k'): xc_pbe0, \
         ('pbe','nlcc_aw'): xc_pbe_nlcc_aw,('pbe','nlcc_ss'): xc_pbe_nlcc_ss}


# Routines for data analysis

def eval_relative_error(dataset,molecule,xc,psp):
    """"
    Compute the relative error associated to the selected molecule in given
    study. Return an array built as
    eps_i = (alpha_ii -alpha_ref_ii)/alpha_ref_ii * 100
    If the reference results for the selected xc and/or the the results for
    the given study (xc+psp) are absent it returns None.
    """
    eps = None
    ref_data = dataset[molecule]['ref_results']
    data = dataset[molecule]['results']
    if xc in ref_data and (xc,psp) in data :
        alpha_ii = get_study_alpha_diag(data,xc,psp)
        alpha_ref_ii = np.array(ref_data[xc])/AuToA
        eps = 100.0 * (alpha_ii-alpha_ref_ii)/alpha_ref_ii

    return eps

def eval_sre_molecule(dataset,mol,xc,psp):
    """
    Compute the (averaged) squared relative error for a single molecule, associated
    to the study xc+psp.
    """
    sre = None
    eps = eval_relative_error(dataset,mol,xc,psp)
    if not (eps is None):
        sre = sum(eps**2)/3.
    return sre

def eval_re_molecule(dataset,mol,xc,psp):
    """
    Compute the (average) relative error for a single molecule, associated
    to the study xc+psp.
    """
    re = None
    eps = eval_relative_error(dataset,mol,xc,psp)
    if not (eps is None):
        re = sum(eps)/3.
    return re

def eval_rmsre(dataset,xc,psp):
    """
    Compute the RMSRE, as defined in the paper of HG, associated
    to the study xc+psp.
    """
    rmsre = 0.
    N = 0
    for mol in build_computed_mol_list(dataset):
        sre = eval_sre_molecule(dataset,mol,xc,psp)
        if not (sre is None):
            rmsre+=sre
            N+=1
    if N>0 : rmsre = rmsre/N
    return np.sqrt(rmsre)

def eval_mre(dataset,xc,psp):
    """
    Compute the MRE, as defined in the paper of HG, associated
    to the study xc+psp.
    """
    mre = 0.
    N = 0
    for mol in build_computed_mol_list(dataset):
        re = eval_re_molecule(dataset,mol,xc,psp)
        if not (re is None):
            mre+=re
            N+=1
    if N>0 : mre = mre/N
    return mre
