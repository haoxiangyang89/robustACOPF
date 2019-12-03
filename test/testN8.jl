# test N8: test the lower bound of infeasibility measure
addprocs(20);
@everywhere using JuMP, MathProgBase, Ipopt, Gurobi, Mosek, Combinatorics;
@everywhere using JLD,HDF5;
@everywhere include("../src/loadMod.jl");

fileAddSet = ["../data/nesta_case5_pjm.m","../data/nesta_case9_wscc.m","../data/nesta_case14_ieee.m","../data/nesta_case30_ieee.m",
    "../data/nesta_case118_ieee.m","../data/nesta_case300_ieee.m","../data/nesta_case2383wp_mp.m","../data/nesta_case2746wp_mp.m"];
nBus = [5,9,14,30,118,300,2383,2746];
ΓSet = [1,3,5];
β = 0.05;
αmax = [0.129,0.213,0.128,0.097,0.140,0.061,0.013,0.033];
dDict = Dict();
groupData = load("../data/groupDict.jld");
groupDict = groupData["groupDict"];
solData = load("../data/testN2Results.jld");
solData = solData["dDict"];

# precompile the function
fi = 1;
fData = makeFData(fileAddSet[fi]);
fData,uData = makeUAData(fData,0,0);
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

for fi in 1:length(fileAddSet)
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
  dDict[fi] = Dict();
  for Γ in ΓSet
    Ω,ypo,yqo = genPartialG(fData,uData,groupDict[fi][2],Γ);
    sphatb,sqhatb = solData[fi][Γ][2];
    sphat0,sqhat0 = solData[fi][Γ][8];
    sphatnc,sqhatnc = solData[fi][Γ][15];

    vioCList0N = testSolC(fData,uData,sphat0,sqhat0,Ω,ypo,yqo,"L",vmax,vmin,θDmax,θDmin)/totalD;
    vioCList0 = [];
    for item in vioCList0N
      if item >= 0
        push!(vioCList0,item);
      else
        push!(vioCList0,0);
      end
    end
    vioCount0 = sum(vioCList0[i] > 1e-4 for i in Ω);
    vioMean0 = mean(vioCList0);
    vioMax0 = maximum(vioCList0);
    dDict[fi][Γ] = [vioCList0,vioCount0,vioMean0,vioMax0];
    println("------------------------- $fi,$Γ finished --------------------------");
    save("../data/testN8Results.jld","dDict",dDict);
  end
end
