# translate the generation at each generator to the injection at each bus
function busInj(fData,sphat,sqhat)
  sphatsum = Dict();
  for i in fData.IDList
    sphatsum[i] = 0.0;
    if i in keys(fData.LocRev)
      for j in fData.LocRev[i]
        sphatsum[i] += sphat[j];
      end
    end
  end

  sqhatsum = Dict();
  for i in fData.IDList
    sqhatsum[i] = 0.0;
    if i in keys(fData.LocRev)
      for j in fData.LocRev[i]
        sqhatsum[i] += sqhat[j];
      end
    end
  end
  return sphatsum,sqhatsum;
end

# generate ∞-norm uncertainty set vertices
function enumN(n,I)
  nTemp = n;
  m = length(I);
  outputList = zeros(m*2);
  counter = 0;
  while nTemp > 1
    outputList[2*m-counter] = mod(nTemp,2);
    nTemp = div(nTemp,2);
    counter += 1;
  end
  if n > 0
    outputList[2*m-counter] = 1;
  else
    outputList[2*m-counter] = 0;
  end
  yps = Dict();
  yqs = Dict();
  for i in 1:length(I)
    yps[i] = outputList[i];
    yqs[i] = outputList[i + length(I)];
  end
  return yps,yqs;
end

# generate ∞-norm uncertainty set vertices
function genOmegaInf(fData)
  Ω = 1:2^(2*length(fData.IDList));
  yp = Dict();
  yq = Dict();
  for ω in Ω
    yp[ω],yq[ω] = enumN(ω - 1, fData.IDList);
  end
  return Ω,yp,yq;
end

# generate 1-norm uncertainty set vertices
function genOmega1(fData)
  # generate the vertices of the 1 norm ball
  Ω = 1:(4*length(fData.IDList));
  yp = 0.5*ones(4*length(fData.IDList),length(fData.IDList));
  yq = 0.5*ones(4*length(fData.IDList),length(fData.IDList));
  ω = 1;
  for i in 1:length(fData.IDList)
    yp[ω,i] = 1;
    ω += 1;
  end
  for i in 1:length(fData.IDList)
    yp[ω,i] = 0;
    ω += 1;
  end
  for i in 1:length(fData.IDList)
    yq[ω,i] = 1;
    ω += 1;
  end
  for i in 1:length(fData.IDList)
    yq[ω,i] = 0;
    ω += 1;
  end
  return Ω,yp,yq;
end

function genPartial1(fData,uData,partCorr = [])
  # generate the vertices of the 1 norm ball with a subset of buses uncertain
  Ωlen = 0;
  yp = Array{Float64}(1,length(fData.IDList));
  yq = Array{Float64}(1,length(fData.IDList));
  iInd = 0;
  dupList = [];
  dupDict = Dict();
  for part in partCorr
    for i in 2:length(part)
      push!(dupList,part[i]);
    end
    dupDict[part[1]] = part[2:length(part)];
  end
  for i in fData.IDList
    if !(i in dupList)
      iInd += 1;
      if uData[i].GPmax - uData[i].GPmin > 1e-4
        # append the scenario where u_p takes the lower bound
        yprow = 0.5*ones(1,length(fData.IDList));
        yqrow = 0.5*ones(1,length(fData.IDList));
        if fData.IDList[iInd] in keys(dupDict)
          for ifData in dupDict[fData.IDList[iInd]]
            iRev = findin(fData.IDList,ifData)[1];
            yprow[iRev] = 0;
          end
        else
          yprow[iInd] = 0;
        end
        yp = [yp;yprow];
        yq = [yq;yqrow];
        Ωlen += 1;
        # append the scenario where u_p takes the upper bound
        yprow = 0.5*ones(1,length(fData.IDList));
        yqrow = 0.5*ones(1,length(fData.IDList));
        if fData.IDList[iInd] in keys(dupDict)
          for ifData in dupDict[fData.IDList[iInd]]
            iRev = findin(fData.IDList,ifData)[1];
            yprow[iRev] = 1;
          end
        else
          yprow[iInd] = 1;
        end
        yp = [yp;yprow];
        yq = [yq;yqrow];
        Ωlen += 1;
      end
      if uData[i].GQmax - uData[i].GQmin > 1e-4
        # append the scenario where u_p takes the lower bound
        yprow = 0.5*ones(1,length(fData.IDList));
        yqrow = 0.5*ones(1,length(fData.IDList));
        if fData.IDList[iInd] in keys(dupDict)
          for ifData in dupDict[fData.IDList[iInd]]
            iRev = findin(fData.IDList,ifData)[1];
            yqrow[iRev] = 0;
          end
        else
          yqrow[iInd] = 0;
        end
        yp = [yp;yprow];
        yq = [yq;yqrow];
        Ωlen += 1;
        # append the scenario where u_p takes the upper bound
        yprow = 0.5*ones(1,length(fData.IDList));
        yqrow = 0.5*ones(1,length(fData.IDList));
        if fData.IDList[iInd] in keys(dupDict)
          for ifData in dupDict[fData.IDList[iInd]]
            iRev = findin(fData.IDList,ifData)[1];
            yqrow[iRev] = 1;
          end
        else
          yqrow[iInd] = 1;
        end
        yp = [yp;yprow];
        yq = [yq;yqrow];
        Ωlen += 1;
      end
    end
  end
  yp = yp[2:(Ωlen+1),:];
  yq = yq[2:(Ωlen+1),:];
  return 1:Ωlen,yp,yq;
end

function sampleY(fData,n)
  Ω = 1:n;
  yp = rand(n,length(fData.IDList));
  yq = rand(n,length(fData.IDList));
  return Ω,yp,yq;
end

function sampleYG(fData,n,partCorr)
  Ω = 1:n;
  partRev = Dict();
  for part in 1:length(partCorr)
    for iBus in partCorr[part]
      partRev[iBus] = part;
    end
  end

  ypr = rand(n,length(partCorr));
  yqr = rand(n,length(partCorr));

  yp = zeros(n,length(fData.IDList));
  yq = zeros(n,length(fData.IDList));
  for i in 1:n
    for j in fData.IDList
      yp[i,fData.busInd[j]] = ypr[i,partRev[j]];
      yq[i,fData.busInd[j]] = yqr[i,partRev[j]];
    end
  end
  return Ω,yp,yq;
end

function sampleYGcorr(fData,n,partCorr)
  Ω = 1:n;
  partRev = Dict();
  for part in 1:length(partCorr)
    for iBus in partCorr[part]
      partRev[iBus] = part;
    end
  end

  ypr = rand(n,length(partCorr));

  yp = zeros(n,length(fData.IDList));
  yq = zeros(n,length(fData.IDList));
  for i in 1:n
    for j in fData.IDList
      yp[i,fData.busInd[j]] = ypr[i,partRev[j]];
      yq[i,fData.busInd[j]] = ypr[i,partRev[j]];
    end
  end
  return Ω,yp,yq;
end

function sampleYGcorrGamma(fData,n,partCorr,Γ)
  Ω = 1:n;
  partRev = Dict();
  for part in 1:length(partCorr)
    for iBus in partCorr[part]
      partRev[iBus] = part;
    end
  end

  nTotal = 0;
  yp = zeros(n,length(fData.IDList));
  yq = zeros(n,length(fData.IDList));
  while nTotal < n
      ypr = rand(n,length(partCorr));
      for i in 1:n
        if (sum(abs.(ypr[i,:] - 0.5)) <= Γ*0.5)&(nTotal < n)
          nTotal += 1;
          for j in fData.IDList
            yp[nTotal,fData.busInd[j]] = ypr[i,partRev[j]];
            yq[nTotal,fData.busInd[j]] = ypr[i,partRev[j]];
          end
        end
      end
  end
  return Ω,yp,yq;
end

function genPartialG(fData,uData,partCorr,b)
  # generate Γ-norm uncertainty set extreme points
  partSize = length(partCorr);
  cPart = combinations(1:partSize,b);
  Ωlen = 0;
  yp = Array{Float64}(1,length(fData.IDList));
  yq = Array{Float64}(1,length(fData.IDList));
  for c in cPart
    for i in 0:(2^b-1)
      yprow = 0.5*ones(1,length(fData.IDList));
      yqrow = 0.5*ones(1,length(fData.IDList));
      # get the binary encoding for i
      iStr = bin(i,b);
      for j in 1:length(iStr)
        jDig = Int64(iStr[j]) - 48;
        for iBus in partCorr[c[j]]
          yprow[fData.busInd[iBus]] = jDig;
          yqrow[fData.busInd[iBus]] = jDig;
        end
      end
      yp = [yp;yprow];
      yq = [yq;yqrow];
      Ωlen += 1;
    end
  end
  yp = yp[2:(Ωlen+1),:];
  yq = yq[2:(Ωlen+1),:];
  return 1:Ωlen,yp,yq;
end

function bestRev(fData,yp,yq)
  # transform the yp/yq from -1/1/0 scale to 0/1/0.5 scale
  yphat = 0.5*ones(1,length(fData.IDList));
  yqhat = 0.5*ones(1,length(fData.IDList));
  for i in fData.IDList
    if yp[i] == 1
      yphat[fData.busInd[i]] = 1;
    elseif yp[i] == -1
      yphat[fData.busInd[i]] = 0;
    end
    if yq[i] == 1
      yqhat[fData.busInd[i]] = 1;
    elseif yq[i] == -1
      yqhat[fData.busInd[i]] = 0;
    end
  end
  return yphat,yqhat;
end

function costCal(fData,sphat,sqhat)
  # calculate the cost of the given solution
  costExpr = 0;
  for i in fData.genIDList
    cpGen = fData.cp[i];
    for j in 1:cpGen.n
      costExpr += cpGen.params[j]*(sphat[i]*fData.baseMVA)^(cpGen.n-j);
    end
    if !(fData.cq == Dict())
      cqGen = fData.cq[i];
      for j in 1:cqGen.n
        costExpr += cqGen.params[j]*(sqhat[i]*fData.baseMVA)^(cqGen.n-j);
      end
    end
  end
  return costExpr
end

function costFix(fData)
  # add a small cost for reactive power generation
  if fData.cq == Dict()
    for i in fData.genIDList
      fData.cq[i] = costDataType(2,0.0,0.0,2,[0.01/fData.baseMVA,0.0]);
    end
  end
  for i in fData.genIDList
    zeroBool = true;
    for j in fData.cp[i].params
      if j != 0.0
        zeroBool = false;
      end
    end
    if zeroBool
      fData.cp[i].params[fData.cp[i].n - 1] = 0.01/fData.baseMVA;
    end
  end
  return fData;
end

function dMatGen(fData)
  # generate the K-cluster based on the network topology
  # obtain the distance between each pair of nodes first
  dMat = zeros(length(fData.IDList),length(fData.IDList));
  # construct the connection matrix
  succDict = Dict();
  for i in fData.IDList
    succDict[i] = [];
  end
  for k in fData.brList
    if !(k[2] in succDict[k[1]])
      push!(succDict[k[1]],k[2]);
      push!(succDict[k[2]],k[1]);
    end
  end

  for i in fData.IDList
    exploredNodes = [i];
    deletedNodes = [];
    while exploredNodes != []
      currentNode = exploredNodes[1];
      shift!(exploredNodes);
      push!(deletedNodes,currentNode);
      for j in succDict[currentNode]
        if (!(j in deletedNodes))&(!(j in exploredNodes))
          push!(exploredNodes,j);
          dMat[fData.busInd[i],fData.busInd[j]] = dMat[fData.busInd[i],fData.busInd[currentNode]] + 1;
        end
      end
    end
  end

  return dMat;
end

function facilityGrouping(fData,dMat,K,d)
  # with dMat, formulate the facility location problem
  facM = Model(solver = GurobiSolver());
  @variable(facM,x[i in fData.IDList, j in fData.IDList;dMat[fData.busInd[j],fData.busInd[i]] <= d],Bin);
  @variable(facM,y[i in fData.IDList],Bin);
  @constraint(facM, connConstr[i in fData.IDList], sum(x[i,j] for j in fData.IDList if dMat[fData.busInd[j],fData.busInd[i]] <= d) == 1);
  @constraint(facM, openConstr[i in fData.IDList, j in fData.IDList; dMat[fData.busInd[j],fData.busInd[i]] <= d], x[i,j] <= y[j]);
  @constraint(facM, noLoc, sum(y[i] for i in fData.IDList) == K);
  @objective(facM, Min, sum(sum(x[i,j]*dMat[fData.busInd[i],fData.busInd[j]] for j in fData.IDList if dMat[fData.busInd[i],fData.busInd[j]] <= d) for i in fData.IDList));

  solve(facM);
  yRes = getvalue(facM[:y]);
  xRes = getvalue(facM[:x]);
  yselect = [];
  for i in fData.IDList
    if yRes[i] > 0
      push!(yselect,i);
    end
  end
  partCorr = [];
  for i in yselect
    newPart = [];
    for j in fData.IDList
      if dMat[fData.busInd[j],fData.busInd[i]] <= d
        if xRes[j,i] == 1.0
          push!(newPart,j);
        end
      end
    end
    push!(partCorr,newPart);
  end

  return partCorr;
end

function facilityGrouping_Heu(fData,dMat,K,d)
  # with dMat, formulate the facility location problem
  facM = Model(solver = GurobiSolver());
  @variable(facM,0 <= x[i in fData.IDList, j in fData.IDList;dMat[fData.busInd[j],fData.busInd[i]] <= d] <= 1);
  @variable(facM,0 <= y[i in fData.IDList] <= 1);
  @constraint(facM, connConstr[i in fData.IDList], sum(x[i,j] for j in fData.IDList if dMat[fData.busInd[j],fData.busInd[i]] <= d) == 1);
  @constraint(facM, openConstr[i in fData.IDList, j in fData.IDList; dMat[fData.busInd[j],fData.busInd[i]] <= d], x[i,j] <= y[j]);
  @constraint(facM, noLoc, sum(y[i] for i in fData.IDList) == K);
  @objective(facM, Min, sum(sum(x[i,j]*dMat[fData.busInd[i],fData.busInd[j]] for j in fData.IDList if dMat[fData.busInd[i],fData.busInd[j]] <= d) for i in fData.IDList));

  solve(facM);
  yRes = getvalue(facM[:y]);
  xRes = getvalue(facM[:x]);
  yselect = [];
  for i in fData.IDList
    if yRes[i] > 0
      push!(yselect,i);
    end
  end

  facM1 = Model(solver = GurobiSolver());
  @variable(facM1,x[i in fData.IDList, j in fData.IDList;yselect[j] > 0], Bin);
  @variable(facM1,y[i in fData.IDList;yselect[i] > 0],Bin);
  @constraint(facM1, connConstr[i in fData.IDList], sum(x[i,j] for j in fData.IDList if yselect[j] > 0) == 1);
  @constraint(facM1, openConstr[i in fData.IDList, j in fData.IDList; yselect[j] > 0], x[i,j] <= y[j]);
  @constraint(facM1, noLoc, sum(y[i] for i in fData.IDList if yselect[i] > 0) == K);
  @objective(facM1, Min, sum(sum(x[i,j]*dMat[fData.busInd[i],fData.busInd[j]] for j in fData.IDList if yselect[j] > 0) for i in fData.IDList));

  solve(facM1);
  yRes = getvalue(facM1[:y]);
  xRes = getvalue(facM1[:x]);
  yselect1 = [];
  for i in fData.IDList
    if yRes[i] > 0
      push!(yselect1,i);
    end
  end

  partCorr = [];
  for i in yselect1
    newPart = [];
    for j in fData.IDList
      if yselect[j] > 0
        if xRes[j,i] == 1.0
          push!(newPart,j);
        end
      end
    end
    push!(partCorr,newPart);
  end

  return partCorr;
end
