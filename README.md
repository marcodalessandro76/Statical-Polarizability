# Impact of the usage of psp in the evaluation of the statical polarizability for molecules

The main aim of this analysis is to study the impact of the usage of the pseudopotentials (psp) on the evalutation
of the statical polarizability for a wide class of molecules given by the dataset of the paper of Head-Gordon (HG).
In order to obtain this information we compare the results produced by different psp (keeping fixed the xc potential)
with the benchmark value provided by an all electron computation performed with the same xc. Results of HG are obtained
exactly in this way and so represents a natural benchmark for our tests (the only possible subtle point is given by the fact
that HG uses a localized basis set, but this fact should not introduce a bias in the evaluation of the statical polarizability
if the completeness of the basis is sufficient). As ulterior benchmark we can also use the results of the Norvergian group
(when available), since they perform all electron computations in a wavelet basis.

The analysis is structured in a set of notebooks built in a sort of top-down approach. 
The (hiercarchically ordered) list of notebooks is: 

* __HG Dataset/Benchmark data__ : perform a parsing of the results of HG and build a dictionary saved in hg_data.yaml. 
  The dictionary contains the reference results for the CCSD(T) method and for the DFT computations with xc = lda,pbe,pbe0.
  Also, for each molecule the type of the study (nsp or sp) is reported.
* __Calculations/Construction of the dataset__ : This notebook reads the hg_data.yaml file and build two dictionaries, one for nsp
  and one for the sp dataset extraced by the HG data. For each molecule, the dictionary contains the reference results
  and the 'study' key. The allowed values of the study type are the tuples: (lda_pt,hgh_k),(lda_pw,hgh_k),(pbe,hgh_k),(pbe0,hgh_k). 
  Moreover, for a subset of molecules also two type of nlcc psp are provided (based on the pbe xc), so for these cases there are two further 
  studies : (pbe,nlcc_aw) and (pbe,nlcc_ss).
* __Calculations/Dataset calculator__ : Read the dataset dictionaries built by 'Construction of the dataset' and compute the statical polarizability
  for all the molecules and study type. Save the results in the files (sp)nsp__results.yaml. This files are the starting point for
  subsequent data analysis. The statical polarizability is computed according to the methods (sp)nsp_workflow defined in the Routines.Py.
  The procedures used in these methods are described in the notebooks 'Single study sp/nsp calculator'. 

* __Calculations/Data analysis__ : 
