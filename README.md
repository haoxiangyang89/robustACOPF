# robustACOPF
Codes and data supplemental files for the paper "Robust Optimization for Electricity Generation"

The materials are contained in three folders:
- src: source codes for algorithms, data input/output and all calculations
- test: test scripts to perform experiments calling functions and generate results in the paper
  - testN1.jl: generates the data for Figure 1
  - testN2.jl: generates Table 1
  - testN3.jl: generates Table 2
  - testN4.jl: generates Table 4
  - testN5.jl: generates Table 3/7
  - testN7.jl: test Calafiore-Campi type solutions
  - testN8/N9.jl: test the convex relaxation lower bound measure
- data: data files used for and generated from tests
  - .m files are the power system data used to construct the network. They are in Matpower data format.
  - groupDict.jld is the clustering results generated from Online Supplement A.
  - boundTightened_NNNN.jld contains the tightened bounds for Case NNNN.
  - other .jld files are outputs from the tests.
