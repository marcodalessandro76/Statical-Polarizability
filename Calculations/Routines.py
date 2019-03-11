from BigDFT.Database import Molecules
from BigDFT import Datasets as D, Calculators as C, Inputfiles as I, Logfiles as lf
import StatPol as SP
import os, sys
from shutil import copyfile

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

def xc_lda(inp,molecule):
    inp.set_xc('LDA')

def xc_pbe(inp,molecule):
    inp.set_xc('PBE')

def xc_pbe0(inp,molecule):
    inp.set_xc('PBE0')
    for atom in get_atoms(molecule):
        key = 'psppar.'+atom
        inp[key]={'Pseudopotential XC': 11}

basepath = os.getcwd()
psp_nlcc_path = basepath+'/psp_nlcc'

print 'Routines basepath :', basepath

def xc_pbe_nlcc(inp,molecule):
    inp.set_xc('PBE')
    for atom in get_atoms(molecule):
        src = psp_nlcc_path+'/psppar.'+atom
        dst = basepath+'/Data/'+molecule+'/pbe-nlcc/psppar.'+atom
        copyfile(src,dst)

set_xc = {('lda','hgh-k') : xc_lda, ('pbe','hgh-k') : xc_pbe,\
          ('pbe0','hgh-k'): xc_pbe0, ('pbe','nlcc'): xc_pbe_nlcc}

def nsp_workflow(alpha_conv=1.0e-2,wf_convergence=1.0e-6,hgrids=0.3,rmult_fine=9.0,term_verb=True,**kwargs):
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

    if not os.path.isdir('Data'): os.mkdir('Data')
    os.chdir('Data')
    if not os.path.isdir(molecule): os.mkdir(molecule)
    os.chdir(molecule)
    path=study[0]+'-'+study[1]
    if not os.path.isdir(path): os.mkdir(path)

    #print ''
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
    results['field_conv']=SP.perform_field_convergence(term_verb=term_verb,rt=alpha_conv,\
    run_dir=path,input=inp,runner=code,posinp=posinp,ppf=SP.eval_alpha)
    f=results['field_conv']['converged_value']

    # rmult convergence
    rmult_list=SP.build_rmult_list([gs_conv,rmult_fine])
    if term_verb : print 'Convergence on rmult'
    results['rmult_conv']=SP.perform_rmult_convergence(term_verb=term_verb,\
    rt=alpha_conv,run_dir=path,intensity=f,rmult=rmult_list,input=inp,\
    runner=code,posinp=posinp,ppf=SP.eval_alpha)

    os.chdir('../../')
    return results
