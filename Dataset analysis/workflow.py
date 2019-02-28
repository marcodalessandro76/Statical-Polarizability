from BigDFT import Datasets as D, Calculators as C, Inputfiles as I, Logfiles as lf
from BigDFT.Database import Molecules
import StatPol as SP
import os

def single_study_workflow(alpha_conv=1.0e-2,wf_convergence=1.0e-6,hgrids=0.3,rmult_fine=9.0,omp=2,mpi=4,term_verb=True,**kwargs):
    """
    Perform the complete workflow to compute the statical polarizability
    of a specific study(molecule+xc+psp)

    Args:
        kwargs['molecule']  : the molecule type
        kwargs['xc']        : the xc functional
        kwargs['psp']       : the pseudopotential
    """
    molecule=kwargs['molecule']
    xc = kwargs['xc']
    psp = kwargs['psp']

    study = {}

    if not os.path.isdir(molecule): os.mkdir(molecule)
    os.chdir(molecule)
    path=xc+'-'+psp
    if not os.path.isdir(path): os.mkdir(path)

    print ''
    print 'Compute alpha for : ', molecule, xc, psp

    posinp=Molecules.Molecule(molecule)

    gs_rtol=10*wf_convergence #rel tol for the gs convergence (using the total energy as control quantity)

    inp = I.Inputfile()
    inp.set_hgrid(hgrids)
    inp.set_xc(xc.upper())
    inp.set_wavefunction_convergence(wf_convergence)

    #gs convergence
    rmult_coarse = [1.0*i for i in range(3,12)]
    data = []
    code=C.SystemCalculator(omp=omp,mpi_run='mpirun -np '+str(mpi),skip=True,verbose=False)
    for r in rmult_coarse:
        gs_study = D.Dataset(label=molecule+'_GS',run_dir=path,posinp=posinp)
        gs_study.set_postprocessing_function(SP.get_energy)
        inp.set_rmult(coarse=r,fine=rmult_fine)
        idd={'rmult':r}
        gs_study.append_run(id=idd,runner=code,input=inp)
        data.append(gs_study)
    if term_verb : print 'Seek for gs convergence'
    study['gs_conv'] = SP.seek_convergence(rt=gs_rtol,term_verb=term_verb,label='rmult',values=rmult_coarse,data=data)

    if term_verb : print 'Seek for alpha convergence'
    # alpha field intensity convergence
    conv_val = study['gs_conv']['converged_value']
    gslog = 'log-'+data[rmult_coarse.index(conv_val)].names[0]+'.yaml'
    gs = lf.Logfile(path+os.sep+gslog)
    inp.set_rmult(gs.log['dft']['rmult'])
    if term_verb : 'Convergence on the field intensity'
    study['field_conv']=SP.perform_field_convergence(term_verb=term_verb,rt=alpha_conv,run_dir=path,input=inp,runner=code,posinp=posinp,ppf=SP.eval_alpha)
    f=study['field_conv']['converged_value']
    # alpha rmult convergence
    rmult_list=SP.build_rmult_list(gs)
    if term_verb : 'Convergence on rmult'
    study['rmult_conv']=SP.perform_rmult_convergence(term_verb=term_verb,rt=alpha_conv,run_dir=path,intensity=f,rmult=rmult_list,input=inp,runner=code,posinp=posinp,ppf=SP.eval_alpha)

    os.chdir('../')
    return study
