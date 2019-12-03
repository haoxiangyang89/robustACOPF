@everywhere using JuMP, MathProgBase, Ipopt, Gurobi, Mosek;
@everywhere using JLD,HDF5;
#@everywhere include("/home/haoxiang/ccsi-robustopt/ccsi-robustopt/src/0.60Version/loadMod.jl");
@everywhere include("loadMod.jl");

fileAddSet = ["nesta_case5_pjm.m","nesta_case9_wscc.m","nesta_case14_ieee.m","nesta_case30_ieee.m",
              "nesta_case118_ieee.m","nesta_case300_ieee.m","nesta_case2383wp_mp.m","nesta_case2746wp_mp.m"];
nBus = [5,9,14,30,118,300,2383,2746];
groupDict = Dict();

for fi in 1:length(fileAddSet)
    # generate the group and output the results
    fileAdd = fileAddSet[fi];
    fData = makeFData(fileAdd);
    d = 10;
    K = 5;
    dMat = dMatGen(fData);
    partCorr = facilityGrouping(fData,dMat,K,d);
    groupDict[fi] = [dMat,partCorr];
end
save("groupDict.jld","groupDict",groupDict);
