# this is the function to execute the solving process - parallel multi-cut version
# the output will be the optimal solution

function solveProcess_OnePMR(fData::fixedData,uData::Dict{Any,Any},Ω,yps,yqs,vmax,vmin,θDmax,θDmin,ρ,ϵ = 1e-4,cutNo = 10,cMax = 500)
  # create the initial master/sub problem
  # without any scenario and cuts yet
  mp = createM(fData,uData);
  totalD = sum(abs(fData.Pd[j]) + abs(fData.Qd[j]) for j in fData.IDList);
  # number of cuts so far and the set of cut coefficients
  nC = 0;
  C = [];

  solve(mp);
  lbCost = getobjectivevalue(mp);
  spbest = Dict();
  sqbest = Dict();

  # record the value of decision variable sp and sq
  sphat = Dict();
  sqhat = Dict();
  for i in fData.genIDList
    sphat[i] = getvalue(mp[:sp][i]);
    sqhat[i] = getvalue(mp[:sq][i]);
  end
  ubCost = costCal(fData,sphat,sqhat);

  # obtain the objective value to test whether the current solution is feasible
  maxVio = -99999;
  yp = SharedArray(yps);
  yq = SharedArray(yqs);
  spInfoList = paraSolve_Multiple(fData,uData,Ω,yps,yqs,sphat,sqhat,vmax,vmin,θDmax,θDmin);
  # sort the spInfoList first
  sort!(spInfoList, by = x -> x[3], rev = true);
  if spInfoList[1][3] > maxVio
    maxVio = spInfoList[1][3];
  end
  nonstop = true;

  while (nonstop)&(nC <= cMax)
    while (maxVio > ϵ*totalD)&(nC <= cMax)
      nC += 1;
      println("------------------------------------------------------------------------------------------------------")
      println("Iteration No. = $(nC), violation = $(maxVio)");
      println("------------------------------------------------------------------------------------------------------")

      # select the cutNo most violated cuts
      # detect whether there are more than cutNo scenarios with violation >= ϵ
      sCounter = 1;
      while (sCounter < length(Ω))&(spInfoList[sCounter][3] > ϵ*totalD)
        sCounter += 1;
      end
      if spInfoList[sCounter][3] <= ϵ*totalD
        sCounter -= 1;
      end
      # if there is more than cutNo scenarios violating, take the top cutNo
      if sCounter > cutNo
        sCounter = cutNo;
      end

      for sc in 1:sCounter
        # obtain the dual variables to form a cut
        sp = spInfoList[sc];
        λpr = Dict();
        λqr = Dict();
        for i in fData.genIDList
          λpr[i] = sp[1][i];
          λqr[i] = sp[2][i];
        end
        # append the cut to the master program
        cutTemp = @constraint(mp,maxVio - sum(λpr[i]*(mp[:sp][i] - sphat[i]) for i in fData.genIDList) -
          sum(λqr[i]*(mp[:sq][i]-sqhat[i]) for i in fData.genIDList) <= 0);
        push!(C,cutTemp);
      end

      # solve the master program, repeat the iteration
      mpCopy = copy(mp);
      mpobj = mpCopy.obj;
      for i in fData.genIDList
        mpobj += ρ*((mpCopy[:sp][i] - sphat[i])^2 + (mpCopy[:sq][i] - sqhat[i])^2);
      end
      @objective(mpCopy,Min,mpobj);
      solve(mpCopy);
      lbCost = getobjectivevalue(mpCopy);
      # record the value of decision variable sp and sq
      sphat = Dict();
      sqhat = Dict();
      for i in fData.genIDList
        sphat[i] = getvalue(mpCopy[:sp][i]);
        sqhat[i] = getvalue(mpCopy[:sq][i]);
      end
      ubCost = costCal(fData,sphat,sqhat);

      # obtain the objective value to test whether the current solution is feasible
      maxVio = -99999;
      spInfoList = paraSolve_Multiple(fData,uData,Ω,yps,yqs,sphat,sqhat,vmax,vmin,θDmax,θDmin);
      # sort the spInfoList first
      sort!(spInfoList, by = x -> x[3], rev = true);
      if spInfoList[1][3] > maxVio
        maxVio = spInfoList[1][3];
      end
    end

    # compute the lower bound
    solve(mp);
    lbCost = getobjectivevalue(mp);
    if abs(ubCost - getobjectivevalue(mp)) < 1e-4*ubCost
      nonstop = false;
      for i in fData.genIDList
        spbest[i] = sphat[i];
        sqbest[i] = sqhat[i];
      end
    else
      nonstop = true;
      # record the value of decision variable sp and sq
      sphat = Dict();
      sqhat = Dict();
      for i in fData.genIDList
        sphat[i] = getvalue(mp[:sp][i]);
        sqhat[i] = getvalue(mp[:sq][i]);
      end

      # obtain the objective value to test whether the current solution is feasible
      maxVio = -99999;
      spInfoList = paraSolve_Multiple(fData,uData,Ω,yps,yqs,sphat,sqhat,vmax,vmin,θDmax,θDmin);
      # sort the spInfoList first
      sort!(spInfoList, by = x -> x[3], rev = true);
      if spInfoList[1][3] > maxVio
        maxVio = spInfoList[1][3];
      end
    end
  end

  # sphat/sqhat is the final solution
  return spbest,sqbest,nC,C,maxVio,lbCost,ubCost;
end

function solveProcess_OnePM_New(fData::fixedData,uData::Dict{Any,Any},Ω,yps,yqs,vmax,vmin,θDmax,θDmin,nT = 10,ϵ = 1e-4,cutNo = 10,cMax = 50,numericOpt = 3)
  # create the initial master/sub problem
  # without any scenario and cuts yet
  mp = createM(fData,uData);
  totalD = sum(abs(fData.Pd[j]) + abs(fData.Qd[j]) for j in fData.IDList);
  # number of cuts so far and the set of cut coefficients
  nC = 0;
  C = [];
  yDict = Dict();
  fixedy = [];

  solve(mp);
  lbCost = getobjectivevalue(mp);
  spbest = Dict();
  sqbest = Dict();

  # record the value of decision variable sp and sq
  sphat = Dict();
  sqhat = Dict();
  for i in fData.genIDList
    sphat[i] = getvalue(mp[:sp][i]);
    sqhat[i] = getvalue(mp[:sq][i]);
  end

  # obtain the objective value to test whether the current solution is feasible
  maxVio = -99999;
  yp = SharedArray(yps);
  yq = SharedArray(yqs);
  spInfoList = paraSolve_Multiple(fData,uData,Ω,yp,yq,sphat,sqhat,vmax,vmin,θDmax,θDmin,3);
  # sort the spInfoList first
  sort!(spInfoList, by = x -> x[3], rev = true);
  if spInfoList[1][3] > maxVio
    maxVio = spInfoList[1][3];
  end

  while (maxVio > ϵ*totalD)&(nC <= cMax)
    nC += 1;
    println("------------------------------------------------------------------------------------------------------")
    println("Iteration No. = $(nC), violation = $(maxVio)");
    println("------------------------------------------------------------------------------------------------------")

    # select the cutNo most violated cuts
    # detect whether there are more than cutNo scenarios with violation >= ϵ
    sCounter = 1;
    while (sCounter < length(Ω))&(spInfoList[sCounter][3] > ϵ*totalD)
      sCounter += 1;
    end
    if spInfoList[sCounter][3] <= ϵ*totalD
      sCounter -= 1;
    end
    # if there is more than cutNo scenarios violating, take the top cutNo
    if sCounter > cutNo
      sCounter = cutNo;
    end

    for sc in 1:sCounter
      # obtain the dual variables to form a cut
      sp = spInfoList[sc];
      λpr = Dict();
      λqr = Dict();
      for i in fData.genIDList
        λpr[i] = sp[1][i];
        λqr[i] = sp[2][i];
      end
      # append the cut to the master program
      cutTemp = @constraint(mp,sp[3] - sum(λpr[i]*(mp[:sp][i] - sphat[i]) for i in fData.genIDList) -
        sum(λqr[i]*(mp[:sq][i]-sqhat[i]) for i in fData.genIDList) <= 0);
      push!(C,(sp[3],λpr,λqr,sphat,sqhat,(yps[sp[4],:],yqs[sp[4],:])));
      if sp[4] in keys(yDict)
        yDict[sp[4]] += 1;
      else
        yDict[sp[4]] = 1;
      end
    end

    # solve the master program, repeat the iteration
    mpChanged = false;
    for iKey in keys(yDict)
      if (yDict[iKey] > nT)&(!(iKey in fixedy))
        push!(fixedy,iKey);
        mpChanged = true;
      end
    end
    if mpChanged
      mp = appendScen(fData,uData,1:length(fixedy),yps[fixedy,:],yqs[fixedy,:],vmax,vmin,θDmax,θDmin,C);
      Ωtemp = [];
      for ω in Ω
        if !(ω in fixedy)
          push!(Ωtemp,ω);
        end
      end
      Ω = copy(Ωtemp);
    end
    solve(mp);
    lbCost = getobjectivevalue(mp);
    # record the value of decision variable sp and sq
    sphat = Dict();
    sqhat = Dict();
    for i in fData.genIDList
      sphat[i] = getvalue(mp[:sp][i]);
      sqhat[i] = getvalue(mp[:sq][i]);
    end

    if Ω != []
      # obtain the objective value to test whether the current solution is feasible
      maxVio = -99999;
      spInfoList = paraSolve_Multiple(fData,uData,Ω,yp,yq,sphat,sqhat,vmax,vmin,θDmax,θDmin,3);
      # sort the spInfoList first
      sort!(spInfoList, by = x -> x[3], rev = true);
      if spInfoList[1][3] > maxVio
        maxVio = spInfoList[1][3];
      end
    else
      maxVio = 0;
    end
  end

  # compute the lower bound
  solve(mp);
  lbCost = getobjectivevalue(mp);
  spbest = Dict();
  sqbest = Dict();
  for i in fData.genIDList
    spbest[i] = getvalue(mp[:sp][i]);
    sqbest[i] = getvalue(mp[:sq][i]);
  end

  # sphat/sqhat is the final solution
  return spbest,sqbest,nC,C,maxVio,lbCost;
end
