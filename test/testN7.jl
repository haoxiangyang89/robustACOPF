# test N7: solve the Calafiore/Campi and then test them at the extreme points
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
  Ω,ypo,yqo = genPartialG(fData,uData,groupDict[fi][2],5);
  N = 1000;

  # simulate 1000 samples to solve the Calafiore-Campi
  countList = Array{Int64}(0);
  meanList = Array{Float64}(0);
  maxList = Array{Float64}(0);
  countListC = Array{Int64}(0);
  meanListC = Array{Float64}(0);
  maxListC = Array{Float64}(0);

  ncList = Array{Float64}(0);
  spList = [];
  sqList = [];
  vioNCList = [];
  nNCList = [];
  for nn in 1:20
    Ωt,ypt,yqt = sampleYGcorr(fData,N,groupDict[fi][2]);
    spccnc,sqccnc,ncnccC,vioNC,nNC = solveENC_Proc(fData,uData,Ωt,ypt,yqt,vmax,vmin,θDmax,θDmin,1e-6);
    push!(ncList,ncnccC);
    push!(spList,spccnc);
    push!(sqList,sqccnc);
    push!(vioNCList,vioNC);
    push!(nNCList,nNC); # to obtain the S_N

    vioNCListcc = testSolNC(fData,uData,spccnc,sqccnc,Ω,ypo,yqo,"L",vmax,vmin,θDmax,θDmin)/totalD;
    vioNCListcc = [max(item,0) for item in vioNCListcc];
    vioCountcc = sum(vioNCListcc[i] > 1e-4 for i in Ω);
    vioMeancc = mean(vioNCListcc);
    vioMaxcc = maximum(vioNCListcc);
    push!(countList,vioCountcc);
    push!(meanList,vioMeancc);
    push!(maxList,vioMaxcc);

    vioCListcc = testSolC(fData,uData,spccnc,sqccnc,Ω,ypo,yqo,"L",vmax,vmin,θDmax,θDmin)/totalD;
    vioCListcc = [max(item,0) for item in vioCListcc];
    vioCCountcc = sum(vioCListcc[i] > 1e-4 for i in Ω);
    vioCMeancc = mean(vioCListcc);
    vioCMaxcc = maximum(vioCListcc);
    push!(countListC,vioCCountcc);
    push!(meanListC,vioCMeancc);
    push!(maxListC,vioCMaxcc);
  end
  dDict[fi] = [ncList,countList,meanList,maxList,countListC,meanListC,maxListC,vioNCList,nNCList];
  println("------------------------- $fi finished --------------------------");
  save("../data/testN7Results.jld","dDict",dDict);
end
