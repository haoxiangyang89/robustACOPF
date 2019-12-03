# test N5: for Case 118/300, compare the effectiveness of different ρ/nT
addprocs(20);
@everywhere using JuMP, MathProgBase, Ipopt, Gurobi, Mosek, Combinatorics;
@everywhere using JLD,HDF5;
@everywhere include("../src/loadMod.jl");

Γ = 3;
β = 0.05;
ωp = β;
ωq = β;
αgtop = 1;
αgbot = 1;
penPer = 0.05;
γ = 0.98;
ρList = [0,0.1,1,10];
nTList = [1,2,3,4,5];

groupData = load("../data/groupDict.jld");
groupDict = groupData["groupDict"];

fileAddSet = ["../data/nesta_case5_pjm.m","../data/nesta_case9_wscc.m","../data/nesta_case14_ieee.m","../data/nesta_case30_ieee.m",
    "../data/nesta_case118_ieee.m","../data/nesta_case300_ieee.m","../data/nesta_case2383wp_mp.m","../data/nesta_case2746wp_mp.m"];
nBus = [5,9,14,30,118,300,2383,2746];
β = 0.05;
αmax = [0.129,0.213,0.128,0.097,0.140,0.061,0.013,0.033];

# prerun the tests
fi = 1;
fData = makeFData(fileAddSet[fi]);
fData,uData = makeUAData(fData,0,0);
Ω,ypo,yqo = genPartialG(fData,uData,groupDict[fi][2],1);
vmax = Dict();
vmin = Dict();
θDmax = Dict();
θDmin = Dict();
for i in fData.IDList
  vmax[i] = fData.Vmax[i];
  vmin[i] = fData.Vmin[i];
end
for k in fData.brList
  θDmax[k] = pi/6;
  θDmin[k] = -pi/6;
end
sphatb,sqhatb,nC,C,vio,cbC = solveProcess_Budget_CorrB(fData,uData,vmax,vmin,θDmax,θDmin,groupDict[fi][2],1,1e-4,5,0);
spbest,sqbest,nC,vioSPList,lbList,ubList,tList = solveProcess_BudgetR_test(fData,uData,vmax,vmin,θDmax,θDmin,0,groupDict[1][2],1,1e-4,500,3);

###########################################################################################
for fi in 5:6
  # generate fData and uData
  fileAdd = fileAddSet[fi];
  fData = makeFData(fileAdd);
  totalD = sum(abs(fData.Pd[j]) + abs(fData.Qd[j]) for j in fData.IDList);
  dDict[fi] = Dict();

  # load the bounds
  vmax0 = Dict();
  vmin0 = Dict();
  θDmax0 = Dict();
  θDmin0 = Dict();
  for i in fData.IDList
    vmax0[i] = fData.Vmax[i];
    vmin0[i] = fData.Vmin[i];
  end
  for k in fData.brList
    θDmax0[k] = pi/6;
    θDmin0[k] = -pi/6;
  end

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

  # obtain the solution for αmax
  ωp = β;
  ωq = β;
  αtop = αmax[fi];
  αbot = 0.2*αtop;
  αgtop = 1;
  αgbot = 1;
  penPer = 0.05;
  γ = 0.98;
  βunc = readGenCorr(fData,γ,groupDict[fi][2],penPer,αgtop,αgbot);
  fData,uData = makeUAData(fData,ωp,ωq,αtop,αbot,βunc);
  Γ = 3;

  dDictReg = Dict();
  for ρ in ρList
  # for each ρ, run the regularized test
    tic();
    spbest,sqbest,nC,vioSPList,lbList,ubList,tList = solveProcess_BudgetR_test(fData,uData,vmax,vmin,θDmax,θDmin,ρ,groupDict[fi][2],Γ,1e-4,300,0);
    t = toc();
    dDictReg[ρ] = (spbest,sqbest,nC,vioSPList,lbList,ubList,tList,t);
  end

  dDictNew = Dict();
  for nT in nTList
    # for each nT, run the new test
    tic();
    spbest,sqbest,nC,C,vioSPList,lbList,ubList,tList = solveProcess_Budget_Newtest(fData,uData,vmax,vmin,θDmax,θDmin,groupDict[fi][2],Γ,1e-4,nT,0);
    t = toc();
    dDictNew[nT] = (spbest,sqbest,nC,C,vioSPList,lbList,ubList,tList,t);
  end
  save("../data/testN5Results_$(fi).jld","dDictReg",dDictReg,"dDictNew",dDictNew);
end
