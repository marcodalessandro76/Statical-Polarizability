from BigDFT.Database import Molecules
from BigDFT import Datasets as D, Calculators as C, Inputfiles as I, Logfiles as lf
import StatPol as SP
import os, sys, numpy as np
from shutil import copyfile

AuToA = 0.5291772085**3

def get_atoms(mol):
    """
    Return a list of string with all the atoms of a molecule
    """
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
    # remove duplicates (and loose ordering)
    atoms = list(set(atoms))
    return atoms

def molecule_inlist(mol,psp_list):
    """
    Check if the atoms of a molecule belong to the list psp_list
    Called to check if the nlcc psp can be computed for the molecule
    """
    atoms = get_atoms(mol)
    inlist = True
    for a in atoms:
        if a not in psp_list:
            inlist = False
            break
    return inlist

def xc_lda_pt(inp,molecule):
    inp.set_xc(1)

def xc_lda_pw(inp,molecule):
    inp.set_xc('-001012')
    for atom in get_atoms(molecule):
        key = 'psppar.'+atom
        inp[key]={'Pseudopotential XC': 1}

def xc_pbe(inp,molecule):
    inp.set_xc(11)

def xc_pbe0(inp,molecule):
    inp.set_xc('PBE0')
    for atom in get_atoms(molecule):
        key = 'psppar.'+atom
        inp[key]={'Pseudopotential XC': 11}

basepath = os.getcwd()
psp_nlcc_aw_path = basepath+'/psp_nlcc_aw'
psp_nlcc_ss_path = basepath+'/psp_nlcc_ss'

#print basepath

def xc_pbe_nlcc_aw(inp,molecule):
    inp.set_xc(11)
    for atom in get_atoms(molecule):
        src = psp_nlcc_aw_path+'/psppar.'+atom
        dst = basepath+'/Data/'+molecule+'/pbe-nlcc_aw/psppar.'+atom
        #dst = basepath+'/Data/'+molecule+'/pbe-nlcc_aw-test/psppar.'+atom
        copyfile(src,dst)

def xc_pbe_nlcc_ss(inp,molecule):
    inp.set_xc(11)
    for atom in get_atoms(molecule):
        src = psp_nlcc_ss_path+'/psppar.'+atom
        dst = basepath+'/Data/'+molecule+'/pbe-nlcc_ss/psppar.'+atom
        #dst = basepath+'/Data/'+molecule+'/pbe-nlcc_ss-test/psppar.'+atom
        copyfile(src,dst)

set_xc = {('lda_pt','hgh_k') : xc_lda_pt,('lda_pw','hgh_k') : xc_lda_pw, \
         ('pbe','hgh_k') : xc_pbe, ('pbe0','hgh_k'): xc_pbe0, \
         ('pbe','nlcc_aw'): xc_pbe_nlcc_aw,('pbe','nlcc_ss'): xc_pbe_nlcc_ss}

def nsp_workflow(alpha_conv=1.0e-2,wf_convergence=1.0e-6,hgrids=0.3,\
                 rmult_fine=9.0,term_verb=True,data_folder='Data',**kwargs):
    """
    Perform the complete nsp workflow to compute the statical polarizability
    of a molecule in a specific study(xc+psp)

    Args:
        kwargs['molecule']  : the molecule type
        kwargs['study']     : the tuple (xc,psp)
        kwargs['code']      : the instance of SystemCalculator
    """
    molecule=kwargs['molecule']
    study = kwargs['study']
    code = kwargs['code']
    code.update_global_options(verbose=term_verb)

    results = {}

    if not os.path.isdir(data_folder): os.mkdir(data_folder)
    os.chdir(data_folder)
    if not os.path.isdir(molecule): os.mkdir(molecule)
    os.chdir(molecule)
    path=study[0]+'-'+study[1]
    if not os.path.isdir(path): os.mkdir(path)

    print 'Compute alpha for : ', molecule,study[0], study[1]

    posinp=Molecules.Molecule(molecule)

    #rel tol for the gs convergence (using the total energy as control quantity)
    gs_rtol=10*wf_convergence

    inp = I.Inputfile()
    inp.set_hgrid(hgrids)
    set_xc[study](inp,molecule)
    inp.set_wavefunction_convergence(wf_convergence)

    if term_verb : print 'Seek for gs convergence'

    rmult_coarse = [1.0*i for i in range(3,12)]
    data = []
    for r in rmult_coarse:
        gs_study = D.Dataset(label=molecule+'_GS',run_dir=path,posinp=posinp)
        gs_study.set_postprocessing_function(SP.get_energy)
        inp.set_rmult(coarse=r,fine=rmult_fine)
        idd={'rmult':r}
        gs_study.append_run(id=idd,runner=code,input=inp)
        data.append(gs_study)

    results['gs_conv'] = SP.seek_convergence(rt=gs_rtol,term_verb=term_verb,\
                       label='rmult',values=rmult_coarse,data=data)

    if term_verb : print 'Seek for alpha convergence'

    # field intensity convergence
    gs_conv = results['gs_conv']['converged_value']
    inp.set_rmult(coarse=gs_conv,fine=rmult_fine)
    if term_verb : print 'Convergence on the field intensity'
    results['field_conv']=SP.perform_field_convergence(term_verb=term_verb,field_int=[1e-2,5e-3,1e-3,5e-4,1e-4],\
                         rt=alpha_conv,run_dir=path,input=inp,runner=code,posinp=posinp,ppf=SP.eval_alpha)
    f=results['field_conv']['converged_value']

    # rmult convergence
    rmult_list=SP.build_rmult_list([gs_conv,rmult_fine])
    if term_verb : print 'Convergence on rmult'
    results['rmult_conv']=SP.perform_rmult_convergence(term_verb=term_verb,\
    rt=alpha_conv,run_dir=path,intensity=f,rmult=rmult_list,input=inp,\
    runner=code,posinp=posinp,ppf=SP.eval_alpha)

    os.chdir('../../')
    return results

# Routines for data analysis

def build_computed_mol_list(dataset):
    """
    Build the list of molecules for which the 'results' key is provided
    """
    computed_mol = []
    for mol,data in dataset.iteritems():
        if 'results' in data:
            computed_mol.append(mol)
    computed_mol.sort()
    return computed_mol

def get_study_alpha_diag(data,xc,psp):
    """
    Extract the diagonal part of the polarizability associated to the study xc,psp taking the
    dictionary data as input. If xcx,psp is not in the computed study returns None.

    Args:
        data : is of the form nsp_results[molecule]['results']
    """
    study = data.keys()
    alpha = None
    if (xc,psp) in study:
        results = data[xc,psp]['rmult_conv']
        rmult = results['converged_value']
        alpha = results['results'][rmult]
        alpha = np.diag(alpha)
    return alpha

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
