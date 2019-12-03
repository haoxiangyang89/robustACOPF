function solveProcess_Budget_New(fData::fixedData,uData::Dict{Any,Any},vmax,vmin,θDmax,θDmin,Γ = 1,ϵ = 1e-4,nT = 10,numericOpt = 0,PresolveOpt = -1)
  # create the initial master/sub problem
  # without any scenario and cuts yet
  totalD = sum(abs(fData.Pd[j]) + abs(fData.Qd[j]) for j in fData.IDList);
  mp = createM(fData,uData,numericOpt);
  # initialze an empty list of scenarios to be added to the master problem
  Ω = [];
  noΩ = 0;
  ypo = Array{Float64, 2}[];
  yqo = Array{Float64, 2}[];
  yDict = Dict();
  fixedy = [];

  sp,spobjo = createS_Budget(fData,uData,vmax,vmin,θDmax,θDmin,Γ,fixedy,numericOpt,PresolveOpt);
  # number of cuts so far and the set of cut coefficients
  nC = 0;
  C = [];

  solve(mp);
  lbCost = getobjectivevalue(mp);
  ubCost = Inf;

  # record the value of decision variable sp and sq
  sphat = Dict();
  sqhat = Dict();
  spobj = copy(spobjo);
  for i in fData.genIDList
    sphat[i] = getvalue(mp[:sp][i]);
    sqhat[i] = getvalue(mp[:sq][i]);
    # update the subproblem objective function with the current master solutions
    spobj -= (sphat[i]*sp[:λpi][fData.Loc[i]] + sqhat[i]*sp[:λqi][fData.Loc[i]]);
  end
  @objective(sp,Max,spobj);
#  sp.solver = GurobiSolver(MIPGap = 0.04,MIQCPMethod = 0);
  solve(sp);

  rmListP = [];
  rmListQ = [];
  for i in fData.IDList
    if uData[i].GPmax - uData[i].GPmin < 1e-4
      push!(rmListP,i);
    end
    if uData[i].GQmax - uData[i].GQmin < 1e-4
      push!(rmListQ,i);
    end
  end

  bestyp = Dict();
  bestyq = Dict();
  for i in fData.IDList
    if !(i in rmListP)
      bestyp[i] = round(getvalue(sp[:ypplus][i])) - round(getvalue(sp[:ypminus][i]));
    else
      bestyp[i] = 0;
    end
    if !(i in rmListQ)
      bestyq[i] = round(getvalue(sp[:yqplus][i])) - round(getvalue(sp[:yqminus][i]));
    else
      bestyq[i] = 0;
    end
  end
  # count the current extreme point
  ypext,yqext = bestRev(fData,bestyp,bestyq);
  if (ypext,yqext) in keys(yDict)
    yDict[(ypext,yqext)] += 1;
  else
    yDict[(ypext,yqext)] = 1;
  end

  # obtain the objective value to test whether the current solution is feasible
  vioSP = getobjectivevalue(sp);
  # obtain the dual variables to form a cut
  λpr = Dict();
  λqr = Dict();
  for i in fData.genIDList
    λpr[i] = getvalue(sp[:λpi][fData.Loc[i]]);
    λqr[i] = getvalue(sp[:λqi][fData.Loc[i]]);
  end
  # append the cut to the master program
  nC += 1;
  println("------------------------------------------------------------------------------------------------------")
  println("Iteration No. = $(nC), violation = $(vioSP)");
  println("------------------------------------------------------------------------------------------------------")
  cutTemp = @constraint(mp,vioSP - sum(λpr[i]*(mp[:sp][i] - sphat[i]) for i in fData.genIDList) -
    sum(λqr[i]*(mp[:sq][i]-sqhat[i]) for i in fData.genIDList) <= 0);
  push!(C,(vioSP,λpr,λqr,sphat,sqhat,(ypext,yqext)));

  while (vioSP > ϵ*totalD)&(nC <= 500)
    # solve the master program, repeat the iteration
    # add the regularization term
    mpChanged = false;
    for iKey in keys(yDict)
      if (yDict[iKey] > nT)&(!(iKey in fixedy))
        noΩ += 1;
        push!(Ω,noΩ);
        ypo = [ypo;iKey[1]];
        yqo = [yqo;iKey[2]];
        push!(fixedy,(iKey[1],iKey[2]));
        mpChanged = true;
      end
    end
    if mpChanged
      mp = appendScen(fData,uData,Ω,ypo,yqo,vmax,vmin,θDmax,θDmin,C,numericOpt);
      sp,spobjo = createS_Budget(fData,uData,vmax,vmin,θDmax,θDmin,Γ,fixedy,numericOpt,PresolveOpt);
    end
    solve(mp);
    lbCost = getobjectivevalue(mp);
    # record the value of decision variable sp and sq
    sphat = Dict();
    sqhat = Dict();
    spobj = copy(spobjo);
    for i in fData.genIDList
      sphat[i] = getvalue(mp[:sp][i]);
      sqhat[i] = getvalue(mp[:sq][i]);
      # update the subproblem objective function with the current master solutions
      spobj -= (sphat[i]*sp[:λpi][fData.Loc[i]] + sqhat[i]*sp[:λqi][fData.Loc[i]]);
    end
    @objective(sp,Max,spobj);
    solve(sp);

    rmListP = [];
    rmListQ = [];
    for i in fData.IDList
      if uData[i].GPmax - uData[i].GPmin < 1e-4
        push!(rmListP,i);
      end
      if uData[i].GQmax - uData[i].GQmin < 1e-4
        push!(rmListQ,i);
      end
    end

    bestyp = Dict();
    bestyq = Dict();
    for i in fData.IDList
      if !(i in rmListP)
        bestyp[i] = round(getvalue(sp[:ypplus][i])) - round(getvalue(sp[:ypminus][i]));
      else
        bestyp[i] = 0;
      end
      if !(i in rmListQ)
        bestyq[i] = round(getvalue(sp[:yqplus][i])) - round(getvalue(sp[:yqminus][i]));
      else
        bestyq[i] = 0;
      end
    end

    # count the current extreme point
    ypext,yqext = bestRev(fData,bestyp,bestyq);
    if (ypext,yqext) in keys(yDict)
      yDict[(ypext,yqext)] += 1;
    else
      yDict[(ypext,yqext)] = 1;
    end

    # obtain the objective value to test whether the current solution is feasible
    vioSP = getobjectivevalue(sp);
    # obtain the dual variables to form a cut
    if vioSP > ϵ*totalD
      λpr = Dict();
      λqr = Dict();
      for i in fData.genIDList
        λpr[i] = getvalue(sp[:λpi][fData.Loc[i]]);
        λqr[i] = getvalue(sp[:λqi][fData.Loc[i]]);
      end
      # append the cut to the master program
      nC += 1;
      println("------------------------------------------------------------------------------------------------------")
      println("Iteration No. = $(nC), violation = $(vioSP)");
      println("------------------------------------------------------------------------------------------------------")
      cutTemp = @constraint(mp,vioSP - sum(λpr[i]*(mp[:sp][i] - sphat[i]) for i in fData.genIDList) -
        sum(λqr[i]*(mp[:sq][i]-sqhat[i]) for i in fData.genIDList) <= 0);
      push!(C,(vioSP,λpr,λqr,sphat,sqhat,(ypext,yqext)));
    end
  end

  solve(mp);
  lbCost = getobjectivevalue(mp);
  spbest = Dict();
  sqbest = Dict();
  for i in fData.genIDList
    spbest[i] = getvalue(mp[:sp][i]);
    sqbest[i] = getvalue(mp[:sq][i]);
  end

  # sphat/sqhat is the final solution
  return spbest,sqbest,nC,C,vioSP,lbCost;
end
