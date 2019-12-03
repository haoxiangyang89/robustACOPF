# test N3: for all test cases, simulate 1000 scenarios within the box set, test the
# feasibility violation of the robust convex solution against those scenarios
addprocs(20);
@everywhere using JuMP, MathProgBase, Ipopt, Gurobi, Mosek, Combinatorics;
@everywhere using JLD,HDF5;
@everywhere include("../src/loadMod.jl");

fileAddSet = ["../data/nesta_case5_pjm.m","../data/nesta_case9_wscc.m","../data/nesta_case14_ieee.m","../data/nesta_case30_ieee.m",
    "../data/nesta_case118_ieee.m","../data/nesta_case300_ieee.m","../data/nesta_case2383wp_mp.m","../data/nesta_case2746wp_mp.m"];
nBus = [5,9,14,30,118,300,2383,2746];
β = 0.05;
ωp = β;
ωq = β;
αmax = [0.129,0.213,0.128,0.097,0.140,0.061,0.013,0.033];

N = 1000;
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
  θDmax[k] = pi/6;
  θDmin[k] = -pi/6;
end
# load the solution data
solDict = load("../testN2Results.jld");
solDict = solDict["dDict"];

for fi in 1:length(fileAddSet)
  αtop = αmax[fi];
  αbot = 0.2*αtop;
  αgtop = 1;
  αgbot = 1;
  penPer = 0.05;
  γ = 0.98;
  # generate fData and uData
  fileAdd = fileAddSet[fi];
  #fData = makeFData(join(["/home/haoxiang/ccsi-robustopt/ccsi-robustopt/src/0.60Version/",fileAddSet[fi]]));
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

  # build the uData
  βunc = readGenCorr(fData,γ,groupDict[fi][2],penPer,αgtop,αgbot);
  fData,uData = makeUAData(fData,ωp,ωq,αtop,αbot,βunc);
  Ωt = Dict();
  ypt = Dict();
  yqt = Dict();
  for nn in 1:20
    #Ωt[nn],ypt[nn],yqt[nn] = sampleY(fData,N);
    Ωt[nn],ypt[nn],yqt[nn] = sampleYGcorr(fData,N,groupDict[fi][2]);
  end

  # generate N number of scenarios
  sphatnc = solDict[fi][5][1][1];
  sqhatnc = solDict[fi][5][1][2];
  for Γ in ΓSet
    # read the solution from testN2NResults.jld
    sphatb = solDict[fi][Γ][5][1];
    sqhatb = solDict[fi][Γ][5][2];
    countList = Array{Int64}(0);
    meanList = Array{Float64}(0);
    maxList = Array{Float64}(0);

    countListC = Array{Int64}(0);
    meanListC = Array{Float64}(0);
    maxListC = Array{Float64}(0);

    countListNC = Array{Int64}(0);
    meanListNC = Array{Float64}(0);
    maxListNC = Array{Float64}(0);

    for nn in 1:20
      tic();
      # an alternative here is to sample with the grouping
      vioCListbnc = testSolNC(fData,uData,sphatb,sqhatb,Ωt[nn],ypt[nn],yqt[nn],"L",vmax,vmin,θDmax,θDmin)/totalD;
      vioCListbnc = [max(item,0) for item in vioCListbnc];
      vioCountbnc = 0;
      vioMeanbnc = mean(vioCListbnc);
      vioMaxbnc = maximum(vioCListbnc);
      for ω in Ωt[nn]
        if vioCListbnc[ω] >= 1e-4
          vioCountbnc += 1;
        end
      end
      push!(countList,vioCountbnc);
      push!(meanList,vioMeanbnc);
      push!(maxList,vioMaxbnc);

      vioCListbc = testSolC(fData,uData,sphatb,sqhatb,Ωt[nn],ypt[nn],yqt[nn],"L",vmax,vmin,θDmax,θDmin)/totalD;
      vioCListbc = [max(item,0) for item in vioCListbc];
      vioCountbc = 0;
      vioMeanbc = mean(vioCListbc);
      vioMaxbc = maximum(vioCListbc);
      for ω in Ωt[nn]
        if vioCListbc[ω] >= 1e-4
          vioCountbc += 1;
        end
      end
      push!(countListC,vioCountbc);
      push!(meanListC,vioMeanbc);
      push!(maxListC,vioMaxbc);

      # test the nonconvex robust solution using random samples
      vioNCListnc = testSolNC(fData,uData,sphatnc,sqhatnc,Ωt[nn],ypt[nn],yqt[nn],"L",vmax,vmin,θDmax,θDmin)/totalD;
      vioNCListnc = [max(item,0) for item in vioNCListnc];
      vioCountnc = 0;
      vioMeannc = mean(vioNCListnc);
      vioMaxnc = maximum(vioNCListnc);
      for ω in Ωt[nn]
        if vioNCListnc[ω] >= 1e-4
          vioCountnc += 1;
        end
      end
      push!(countListNC,vioCountnc);
      push!(meanListNC,vioMeannc);
      push!(maxListNC,vioMaxnc);
      tIter = toc();
      println("------------------------- Sample $nn finished, time elapsed $tIter seconds --------------------------");
    end
    dDict[fi][Γ] = [mean(countList),mean(meanList),mean(maxList),std(meanList),std(maxList),
      mean(countListC),mean(meanListC),mean(maxListC),std(meanListC),std(maxListC),
      countList,meanList,maxList,countListNC,meanListNC,maxListNC,countListC,meanListC,maxListC];
    println("------------------------- $fi,$Γ finished --------------------------");
    save("../data/testN3Results.jld","dDict",dDict);
  end

  # output
  open("testN3_Results_$(fi).txt","w") do f
    write(f,"αω, vioCountb, vioMeanb, vioMaxb, sigmaMean, sigmaMax, vioCountbc, vioMeanbc, vioMaxbc, sigmaMeanc, sigmaMaxc\n");
    for Γ in ΓSet
        write(f,"$(Γ),$(dDict[fi][Γ][1]),$(dDict[fi][Γ][2]),$(dDict[fi][Γ][3]),$(dDict[fi][Γ][4]),$(dDict[fi][Γ][5]),
          $(dDict[fi][Γ][6]),$(dDict[fi][Γ][7]),$(dDict[fi][Γ][8]),$(dDict[fi][Γ][9]),$(dDict[fi][Γ][10])\n");
    end
    write(f,"\n");
  end
end
