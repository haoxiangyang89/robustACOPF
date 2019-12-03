# test N2: for all test cases, compare the quality of the robust solution and the
#         nominal solution at the extreme points under nonconvex setting
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
sphatc,sqhatc,ccC,vioC = solveEC_Proc(fData,uData,Ω,ypo,yqo,vmax,vmin,θDmax,θDmin);
sphatnc,sqhatnc,ccnC,vioNC = solveENC_Proc(fData,uData,Ω,ypo,yqo,vmax,vmin,θDmax,θDmin);

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
      sphatnc,sqhatnc,ncnbC,vioNC,nNC = solveENC_Proc(fData,uData,Ω,ypo,yqo,vmax,vmin,θDmax,θDmin);
      tnc = toc();

      # solve for the convex robust solution with bound tightened
      tic();
      sphatb,sqhatb,nC,C,vio,cbC = solveProcess_Budget_CorrB(fData,uData,vmax,vmin,θDmax,θDmin,groupDict[fi][2],Γ,1e-4,1,0);
      tb = toc();
      cbC = costCal(fData,sphatb,sqhatb);
      vioNCListbN = testSolNC(fData,uData,sphatb,sqhatb,Ω,ypo,yqo,"L",vmax0,vmin0,θDmax0,θDmin0)/totalD;
      vioNCListb = [];
      for item in vioNCListbN
        if item >= 0
          push!(vioNCListb,item);
        else
          push!(vioNCListb,0);
        end
      end
      vioCountb = sum(vioNCListb[i] > 1e-4 for i in Ω);
      vioMeanb = mean(vioNCListb);
      vioMaxb = maximum(vioNCListb);

      # test for the deterministic solution
      vioNCList0N = testSolNC(fData,uData,sphat0,sqhat0,Ω,ypo,yqo,"L",vmax0,vmin0,θDmax0,θDmin0)/totalD;
      vioNCList0 = [];
      for item in vioNCList0N
        if item >= 0
          push!(vioNCList0,item);
        else
          push!(vioNCList0,0);
        end
      end
      vioCount0 = sum(vioNCList0[i] > 1e-4 for i in Ω);
      vioMean0 = mean(vioNCList0);
      vioMax0 = maximum(vioNCList0);

      vioCList0N = testSolC(fData,uData,sphat0,sqhat0,Ω,ypo,yqo,"L",vmax0,vmin0,θDmax0,θDmin0)/totalD;
      vioCList0 = [];
      for item in vioCList0N
        if item >= 0
          push!(vioCList0,item);
        else
          push!(vioCList0,0);
        end
      end
      vioCCount0 = sum(vioCList0[i] > 1e-4 for i in Ω);
      vioCMean0 = mean(vioCList0);
      vioCMax0 = maximum(vioCList0);

      # record the results
      dDict[fi][Γ] = [(sphatnc,sqhatnc),ncnbC,tnc,
          cbC,(sphatb,sqhatb),vioNCListb,vioCountb,vioMeanb,vioMaxb,tb,
          cbC0,(sphat0,sqhat0),vioNCList0,vioCount0,vioMean0,vioMax0,vioCList0,vioCCount0,vioCMean0,vioCMax0,t0,
          vioNC,nNC,vioC,nC];
      println("------------------------- $fi,$Γ finished --------------------------");
      save("../data/testN2Results.jld","dDict",dDict);
  end

  # output
  open("testN2_Results_$(fi).txt","w") do f
    write(f,"αω, ncnbC, tnc,
      cbC, vioCountb, vioMeanb, vioMaxb, tb,
      cbC0, vioCount0, vioMean0, vioMax0, vioCCount0, vioCMean0, vioCMax0, t0\n");
    for Γ in ΓSet
        write(f,"$(Γ),$(round(dDict[fi][Γ][2],1)),$(round(dDict[fi][Γ][3],4)),
          $(round(dDict[fi][Γ][4],1)),$(dDict[fi][Γ][7]),$(dDict[fi][Γ][8]),$(dDict[fi][Γ][9]),$(round(dDict[fi][Γ][10],4)),
          $(round(dDict[fi][Γ][11],1)),$(dDict[fi][Γ][14]),$(dDict[fi][Γ][15]),$(dDict[fi][Γ][16]),$(dDict[fi][Γ][18]),$(dDict[fi][Γ][19]),$(dDict[fi][Γ][20]),$(round(dDict[fi][Γ][21],4))
          \n");
    end
    write(f,"\n");
  end
end
