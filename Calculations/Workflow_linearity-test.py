from BigDFT import Datasets as D, Calculators as C, Inputfiles as I, Logfiles as lf
from BigDFT.Database import Molecules
import numpy as np
import sys,os
sys.path.insert(0,'../')
import StatPol as SP

# Service functions
def change_workdir(molecule,study,datadir):
    """
    Ensure that the path Data_directory/molecule exists and
    move there.
    """
    from futile.Utils import ensure_dir
    for d in [datadir,molecule]:
        ensure_dir(d)
        os.chdir(d)
    path='-'.join(study)
    ensure_dir(path)
    return '-'.join([molecule,path]),path

def get_psp(molecule):
    """
    Check if the all the atoms of molecule have the nlcc of the type nlss_aw and nlcc_ss.
    If it is true the corresponding pseudo is added to the list of psp to be performed, otherwise
    only hgh_k is included.
    """
    from os import path as p
    import Routines as R
    possible_psp=['nlcc_aw','nlcc_ss']
    to_do=[True,True]
    for at in R.get_atoms(molecule):
        psp_available = [p.isfile(p.join('psp_'+psp,'psppar.'+at)) for psp in possible_psp]
        to_do = [ a and b for a,b in zip(to_do,psp_available)]
    return ['hgh_k'] + [psp for yes,psp in zip(to_do,possible_psp) if yes]

def split_dataset(dataset,length=10):
    """
    Split a list of molecules into a list of list with length elements.
    """
    splitted = [[] for num in range(int(len(dataset)/length)+1)]
    for ind,mol in enumerate(dataset):
        splitted[int(ind/length)].append(mol)
    return splitted

# Tools for the ground state
def get_converged_input_energy(dataset,rtol,atol):
    data=dataset.seek_convergence(atol=atol,rtol=rtol,attribute='energy')
    return dataset.runs[dataset.names.index(D.name_from_id(data[0]))]['input'],data[1]

def find_gs_domain(label,rtol,atol,path,input,posinp,code):
    """
    Use the seek_convergence method of Dataset to perform a convergence procedure on
    the rmult value for a gs computation.
    """
    crmult=input['dft']['rmult'][0]
    rmult_fine=input['dft']['rmult'][1]
    rmult_list=map(float,range(int(crmult),11))
    seek_for_rmult = D.Dataset(label=label+'(GS)',run_dir=path,posinp=posinp)
    for rm in rmult_list:
        input.set_rmult(coarse=rm,fine=rmult_fine)
        seek_for_rmult.append_run(id={'rmult':rm},runner=code,input=input)
    input_gs,log_gs=get_converged_input_energy(seek_for_rmult,rtol,atol)
    return {'input_gs':input_gs, 'dataset_gs': seek_for_rmult,'log_gs':log_gs}

def gs_study(molecule,study,code,options):
    """"
    Workflow for the convergence analysis of the ground state.

    Args :
       molecule (str) : name of the molecule
       study (touple) : the couple (xc,psp)
       code (runner)  : instance of SystemCalculator
       options (dict) : dictionary with the computational options
    """
    import Routines as R
    hgrids=options.get('hgrids')
    rmult_fine=options.get('rmult_fine')
    wf_convergence=options.get('wf_convergence')
    crmult_start=options.get('crmult',4.0)
    rtol=options.get('rtol_gs',10*wf_convergence)
    atol=options.get('atol_gs',1.e-3)
    reference_results=options.get('reference_data')
    datadir=options.get('data_directory','Data')

    initial_dir=os.path.abspath(os.path.dirname('.'))
    label,path=change_workdir(molecule,study,datadir)

    posinp=Molecules.Molecule(molecule)

    inp = I.Inputfile()
    inp.set_hgrid(hgrids)
    inp.set_rmult(coarse=crmult_start,fine=rmult_fine)
    R.set_xc[study](inp,molecule,datadir)
    inp.set_wavefunction_convergence(wf_convergence)
    mpol=int(reference_results[molecule]['mpol_ref'])-1
    if mpol >0: inp.spin_polarize(mpol)
    data= find_gs_domain(label,10*wf_convergence,atol,path,inp,posinp,code)
    os.chdir(initial_dir)

    data['molecule']=molecule
    data['posinp']=posinp
    data['study']=study
    return data

# Tools for alpha
def extract_statical_polarizability(label,rtol,atol,path,input_for_alpha,posinp,code,ref):

    # seek for the field convergence
    ints=[1e-2,5e-3,1e-3,5e-4,1e-4]
    alpha_field=D.Dataset(label=label+'(field)',run_dir=path)
    for f in ints:
        alpha_field.append_run(id={'f':f},runner=SP.build_alpha_dataset(input=input_for_alpha,posinp=posinp,run_dir=path,runner=code,intensity=f))
    try:
        data_field=alpha_field.seek_convergence(rtol=rtol,atol=atol)
        intensity=data_field[0]['f']
    except:
        print ('Convergence in field not reached')
        intensity=ints[-1]

    # seek for the rmult convergence
    rmult_fine=input_for_alpha['dft']['rmult'][1]
    rmult_list=map(float,range(int(input_for_alpha['dft']['rmult'][0]),11))
    alpha_rmult=D.Dataset(label=label+'(domain)',run_dir=path)
    for rm in rmult_list:
        input_for_alpha.set_rmult(coarse=rm,fine=rmult_fine)
        alpha_rmult.append_run(id={'f':intensity,'rmult':rm},runner=SP.build_alpha_dataset(input=input_for_alpha,posinp=posinp,run_dir=path,runner=code,intensity=intensity))
    try:
        import numpy as np
        data_final=alpha_rmult.seek_convergence(rtol=rtol,atol=atol)
        AuToA = 0.5291772085**3
        alpha_ref_ii = np.array(ref)/AuToA
        eps = 100.0 * (data_final[1]-alpha_ref_ii)/alpha_ref_ii
        print ('Relative difference in %',eps.tolist())
    except LookupError as e:
        print ('Convergence in domain not reached',e)
        data_final=None
    return {'alpha_convergence': data_final, 'dataset_field': alpha_field, 'dataset_rmult': alpha_rmult}

def extract_statical_polarizability_linearity(label,rtol,atol,path,input_for_alpha,posinp,code,ref):

    #ints=[1e-2,5e-3,1e-3,5e-4,1e-4]
    ints= [1e-3,2e-3,3e-3,4e-3,5e-3] #for comparison with the MWL plots

    results = {}
    results['alpha_values'] = {}
    rmult_fine=input_for_alpha['dft']['rmult'][1]
    rmult_list=map(float,range(int(input_for_alpha['dft']['rmult'][0]),12))
    for intensity in ints:

        # seek for the rmult convergence
        alpha_rmult=D.Dataset(label=label+'(domain)',run_dir=path)
        for rm in rmult_list:
            input_for_alpha.set_rmult(coarse=rm,fine=rmult_fine)
            alpha_rmult.append_run(id={'f':intensity,'rmult':rm},\
            runner=SP.build_alpha_dataset(input=input_for_alpha,posinp=posinp,run_dir=path,runner=code,intensity=intensity))

        try:
            import numpy as np
            data_final=alpha_rmult.seek_convergence(rtol=rtol,atol=atol)
            AuToA = 0.5291772085**3
            alpha_ref_ii = np.array(ref)/AuToA
            eps = 100.0 * (data_final[1]-alpha_ref_ii)/alpha_ref_ii
            print ('Relative difference in %',eps.tolist())
        except LookupError as e:
            print ('Convergence in domain not reached',e)
            data_final=None

        results['alpha_values'][intensity] = {'alpha_convergence': data_final, 'dataset_rmult': alpha_rmult}
    return results

def alpha_study(gs_data,code,options):
    """"
    Workflow for the convergence analysis for the computation of alpha.

    Args :
       gs_data        : output of gs_study
       code (runner)  : instance of SystemCalculator
       options (dict) : dictionary with the computational options
    """
    import copy
    xc_conversion={'lda_pw': 'lda-SPW92', 'pbe': 'pbe', 'pbe0': 'pbe0'}
    rtol=options.get('rtol',1.e-2)
    atol=options.get('atol',1.e-3)
    reference_results=options.get('reference_data')
    datadir=options.get('data_directory','Data')
    study=gs_data['study']
    molecule=gs_data['molecule']
    posinp=gs_data['posinp']
    xc=study[0]
    ref_data=reference_results[molecule][xc_conversion[xc]]

    initial_dir=os.path.abspath(os.path.dirname('.'))
    label,path=change_workdir(molecule,study,datadir)
    input_gs=copy.deepcopy(gs_data['input_gs'])

    #data=extract_statical_polarizability(label,rtol,atol,path,input_gs,posinp,code,ref=ref_data)
    data=extract_statical_polarizability_linearity(label,rtol,atol,path,input_gs,posinp,code,ref=ref_data)
    os.chdir(initial_dir)

    return data

def data_to_save(data,options):
    """
    Define the dictionary to be saved on file for data analysis
    """
    from copy import deepcopy
    res = {}
    for s,values in data.iteritems():
        res[s] = {}
        if not (values['alpha_convergence'] is None):
            res[s]['alpha_convergence'] = (values['alpha_convergence'][0],values['alpha_convergence'][1].tolist())
        else :
            res[s]['alpha_convergence'] = values['alpha_convergence']

        res[s]['posinp'] = {}
        pos = values['posinp']
        for key in pos:
            res[s]['posinp'][key] = pos[key]

        res[s]['input_gs'] = {}
        inp = values['input_gs']
        for key in inp:
            res[s]['input_gs'][key] = inp[key]

    save_options = deepcopy(options)
    save_options.pop('reference_data')
    res['options'] = save_options
    return res

def data_to_save_linearity(data,options):
    """
    Define the dictionary to be saved on file for data analysis
    """
    from copy import deepcopy
    res = {}
    for s,values in data.iteritems():
        res[s] = {}
        res[s]['alpha_values'] = {}
        alpha_values = values['alpha_values']
        for int,v in alpha_values.iteritems():
            res[s]['alpha_values'][int] = {}
            if not (v['alpha_convergence'] is None):
                res[s]['alpha_values'][int]['alpha_convergence'] = (v['alpha_convergence'][0],v['alpha_convergence'][1].tolist())
            else:
                res[s]['alpha_values'][int]['alpha_convergence'] = values['alpha_convergence']

        res[s]['posinp'] = {}
        pos = values['posinp']
        for key in pos:
            res[s]['posinp'][key] = pos[key]

        res[s]['input_gs'] = {}
        inp = values['input_gs']
        for key in inp:
            res[s]['input_gs'][key] = inp[key]

    save_options = deepcopy(options)
    save_options.pop('reference_data')
    res['options'] = save_options
    return res

# Multiple study analysis - linearity test

import yaml
HG_data=yaml.load(open('../HG Dataset/hg_data.yaml'))

options={'wf_convergence':1.e-6,'hgrids':0.3,'rmult_fine':9.0,'rtol':1.e-3,'atol':1.e-3,\
'reference_data':HG_data,'data_directory':'Data'}

# mpi and omp are set from above exporting the associated variables
code=C.SystemCalculator(skip=True,verbose=False)

studies = [('Na','pbe','hgh_k'),('Li','pbe','nlcc_ss')]

study_data = {}

for study in studies:
    molecule = study[0]
    xc = study[1]
    psp = study[2]

    gs_data=gs_study(molecule,(xc,psp),code,options)
    study_data[(molecule,xc,psp)]=alpha_study(gs_data,code,options)
    study_data[(molecule,xc,psp)].update(gs_data)

full_results_linearity = {}
#full_results_linearity = yaml.load(open('full_results_linearity.yaml'))

full_results_linearity.update(data_to_save_linearity(study_data,options))
with open('full_results_linearity.yaml', 'w') as outfile:
    yaml.dump(full_results_linearity, outfile, default_flow_style=False)
