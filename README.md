# Impact of the usage of psp in the evaluation of the statical polarizability for molecules

The main aim of this analysis is to study the impact of the usage of the pseudopotentials (psp) on the evalutation
of the statical polarizability for a wide class of molecules given by the dataset of the paper of Head-Gordon (HG).
In order to obtain this information we compare the results produced by different psp (keeping fixed the xc potential)
with the benchmark value provided by an all electron computation performed with the same xc. Results of HG are obtained
exactly in this way and so represents a natural benchmark for our tests (the only possible subtle point is given by the fact
that HG uses a localized basis set, but this fact should not introduce a bias in the evaluation of the statical polarizability
if the completeness of the basis is sufficient). As ulterior benchmark we can also use the results of the Norvergian group
(when available), since they perform all electron computations in a wavelet basis.

Use a sort of top-down approach. Build the ordered structure of notebook:

* Benchmark data
* Construction of the dataset
* Dataset calculator
