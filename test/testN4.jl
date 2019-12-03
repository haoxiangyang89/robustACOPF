# test N4: for all test cases, compare the solution process time of the MISOCP version
# vs. the multiple-SOCP version.

addprocs(20);
@everywhere using JuMP, MathProgBase, Ipopt, Gurobi, Mosek, Combinatorics;
@everywhere using JLD,HDF5;
@everywhere include("../src/loadMod.jl");

fileAddSet = ["../data/nesta_case5_pjm.m","../data/nesta_case9_wscc.m","../data/nesta_case14_ieee.m","../data/nesta_case30_ieee.m",
    "../data/nesta_case118_ieee.m","../data/nesta_case300_ieee.m","../data/nesta_case2383wp_mp.m","../data/nesta_case2746wp_mp.m"];
nBus = [5,9,14,30,118,300];
αmax = [0.129,0.213,0.128,0.097,0.140,0.061,0.013,0.033];
β = 0.05;
αgtop = 1;
αgbot = 1;
penPer = 0.05;
γ = 0.98;
ΓSet = [1,3,5];

dDict = Dict();
groupData = load("../data/groupDict.jld");
groupDict = groupData["groupDict"];

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
  θDmax[k] = pi/3;
  θDmin[k] = -pi/3;
end
sphatb,sqhatb,nC,C,vio,cbC = solveProcess_Budget_CorrB(fData,uData,vmax,vmin,θDmax,θDmin,groupDict[fi][2],1,1e-4,5,3);
# revise the OnePMR to utilize a SharedArray
Ω,ypo,yqo = genPartialG(fData,uData,groupDict[fi][2],1);
sphat2,sqhat2,nC2,C2,vio2,cbC2= solveProcess_OnePM_New(fData,uData,Ω,ypo,yqo,vmax,vmin,θDmax,θDmin,5,1e-4);

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
      θDmax[k] = pi/3;
      θDmin[k] = -pi/3;
    end
  end

  # build the uData
  ωp = β;
  ωq = β;
  αtop = αmax[fi];
  αbot = 0.2*αtop;
  βunc = readGenCorr(fData,γ,groupDict[fi][2],penPer,αgtop,αgbot);
  fData,uData = makeUAData(fData,ωp,ωq,αtop,αbot,βunc);
  dDict[fi] = Dict();
  for Γ in ΓSet
    Ω,ypo,yqo = genPartialG(fData,uData,groupDict[fi][2],Γ);

    # solve the problem using MISOCP
    tic();
    sphatc,sqhatc,nC,C,vio,cbC = solveProcess_Budget_CorrB(fData,uData,vmax,vmin,θDmax,θDmin,groupDict[fi][2],Γ,1e-4,1,3);
    misocpt = toc();

    # solve the problem using SOCP
    tic();
    sphat2,sqhat2,nC2,C2,vio2,cbC2 = solveProcess_OnePM_New(fData,uData,Ω,ypo,yqo,vmax,vmin,θDmax,θDmin,1,1e-4,10,500,3);
    enumt = toc();

    dDict[fi][Γ] = [Γ,length(Ω),misocpt,enumt,nC,nC2,sphatc,sqhatc,sphat2,sqhat2];
    save("../data/testN4Results.jld","dDict",dDict);
  end
  println(" ------------------------- $fi finished -------------------------- ");
end
# output
open("testN4_Results.txt","w") do f
  write(f,"Gamma, Vertices, MISOCP Time, Enumeration Time (Par), MISOCP No. Iteration, SOCP No. Iteration\n");
  for fi in 1:length(fileAddSet)
    for Γ in ΓSet
      write(f,"$(dDict[fi][Γ][1]),$(dDict[fi][Γ][2]),$(dDict[fi][Γ][3]),$(dDict[fi][Γ][4]),$(dDict[fi][Γ][5]),$(dDict[fi][Γ][6])\n");
    end
  end
  write(f,"\n");
end
