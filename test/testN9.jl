# test N9: test the deterministic and extensive convex relaxation
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
  for Γ in ΓSet
      # build the uData
      Ω,ypo,yqo = genPartialG(fData,uData,groupDict[fi][2],Γ);
      tic();
      sphat0,sqhat0,cbC0 = solveNominalC(fData,uData,vmax,vmin,θDmax,θDmin);
      t0 = toc();
      cbC0 = costCal(fData,sphat0,sqhat0);
      tic();
      sphatc,sqhatc,ccC = solveEC_Proc(fData,uData,Ω,ypo,yqo,vmax,vmin,θDmax,θDmin);
      tc = toc();
      ccC = costCal(fData,sphatc,sqhatc);

      tic();
      if fi <= 6
        sphatb,sqhatb,nC,C,vio,cbC = solveProcess_Budget_CorrB(fData,uData,vmax,vmin,θDmax,θDmin,groupDict[fi][2],Γ,1e-4,1,0);
      else
        sphatb,sqhatb,nC,C,vio,cbC = solveProcess_Budget_CorrB(fData,uData,vmax,vmin,θDmax,θDmin,groupDict[fi][2],Γ,1e-4,1,0,0);
      end
      tb = toc();

      vioNCListb = testSolNC(fData,uData,sphatb,sqhatb,Ω,ypo,yqo,"L",vmax,vmin,θDmax,θDmin)/totalD;
      vioNCListb = [max(item,0) for item in vioNCListb];
      vioCountb = sum(vioNCListb[i] > 1e-4 for i in Ω);
      vioMeanb = mean(vioNCListb);
      vioMaxb = maximum(vioNCListb);

      # test for the deterministic solution
      vioNCList0 = testSolNC(fData,uData,sphat0,sqhat0,Ω,ypo,yqo,"L",vmax,vmin,θDmax,θDmin)/totalD;
      vioNCList0 = [max(item,0) for item in vioNCList0];
      vioCount0 = sum(vioNCList0[i] > 1e-4 for i in Ω);
      vioMean0 = mean(vioNCList0);
      vioMax0 = maximum(vioNCList0);

      vioNCListc = testSolNC(fData,uData,sphatc,sqhatc,Ω,ypo,yqo,"L",vmax,vmin,θDmax,θDmin)/totalD;
      vioNCListc = [max(item,0) for item in vioNCListc];
      vioCountc = sum(vioNCListc[i] > 1e-4 for i in Ω);
      vioMeanc = mean(vioNCListc);
      vioMaxc = maximum(vioNCListc);

      # record the results
      dDict[fi][Γ] = [cbC,(sphatb,sqhatb),vioNCListb,vioCountb,vioMeanb,vioMaxb,tb,
          cbC0,(sphat0,sqhat0),vioNCList0,vioCount0,vioMean0,vioMax0,t0,
          ccC,(sphatc,sqhatc),vioNCListc,vioCountc,vioMeanc,vioMaxc,tc];
      println("------------------------- $fi,$Γ finished --------------------------");
      save("../data/testN9Results.jld","dDict",dDict);
  end
end
