# test N1: for small test cases, obtain the alpha_max, solve for the nonconvex robust
#         solution, compare it with the convex robust solution.
addprocs(20);
@everywhere using JuMP, MathProgBase, Ipopt, Gurobi, Mosek, Combinatorics;
@everywhere using JLD,HDF5;
@everywhere include("../src/loadMod.jl");

fileAddSet = ["../data/nesta_case5_pjm.m","../data/nesta_case9_wscc.m","../data/nesta_case14_ieee.m","../data/nesta_case30_ieee.m",
    "../data/nesta_case118_ieee.m","../data/nesta_case300_ieee.m","../data/nesta_case2383wp_mp.m","../data/nesta_case2746wp_mp.m"];
nBus = [5,9,14,30,118,300,2383,2746];
budgetSet = [1,3,5];
βList = [0.01,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4];
dDict = Dict();
groupData = load("../data/groupDict.jld");
groupDict = groupData["groupDict"];

for fi in 1:length(fileAddSet)
  # generate fData and uData
  fileAdd = fileAddSet[fi];
  fData = makeFData(fileAdd);
  totalD = sum(abs(fData.Pd[j]) + abs(fData.Qd[j]) for j in fData.IDList);
  dDict[fi] = Dict();

  # load the bounds
  vmax = Dict();
  vmin = Dict();
  θDmax = Dict();
  θDmin = Dict();
  try
    boundData = load("../data/boundTightened_$(nBus[fi]).jld");
    vmax = boundData["boundDict"][1];
    vmin = boundData["boundDict"][2];
    θDmax = boundData["boundDict"][3];
    θDmin = boundData["boundDict"][4];
  catch
    for i in fData.IDList
      vmax[i] = fData.Vmax[i];
      vmin[i] = fData.Vmin[i];
    end
    for k in fData.brList
      θDmax[k] = pi/6;
      θDmin[k] = -pi/6;
    end
  end

  # obtain αmax for nonconvex, convex, convex nonsymmetric settings
  rDictNC = Dict();
  αmaxListNC = Dict();
  rDictC = Dict();
  αmaxListC = Dict();
  αgtop = 1;
  αgbot = 1;
  penPer = 0.05;
  γ = 0.98;
  βunc = readGenCorr(fData,γ,groupDict[fi][2],penPer,αgtop,αgbot);

  for tmb in budgetSet
    rDictNC[tmb],αmaxListNC[tmb] = getAlphaNC(fData,βList,tmb,vmax,vmin,θDmax,θDmin,groupDict[fi][2],0.2,1e-4,βunc);
    save("../data/alphaMaxResults_NC_$(fi).jld","rDictNC",rDictNC,"αmaxListNC",αmaxListNC);
    rDictC[tmb],αmaxListC[tmb] = getAlphaC(fData,βList,tmb,vmax,vmin,θDmax,θDmin,groupDict[fi][2],0.2,1e-4,βunc);
    save("../data/alphaMaxResults_C_$(fi).jld","rDictC",rDictC,"αmaxListC",αmaxListC);
  end
  println("------------------------- $fi αmax obtained --------------------------");
end
