function solveProcess_BudgetR(fData::fixedData,uData::Dict{Any,Any},vmax,vmin,θDmax,θDmin,ρ,Γ = 1, ϵ = 1e-4, cMax = 500)
  # create the initial master/sub problem
  # without any scenario and cuts yet
  totalD = sum(abs(fData.Pd[j]) + abs(fData.Qd[j]) for j in fData.IDList);
  mp = createM(fData,uData);
  sp,spobjo = createS_Budget(fData,uData,vmax,vmin,θDmax,θDmin,Γ);
  # number of cuts so far and the set of cut coefficients
  nC = 0;
  C = [];

  solve(mp);
  lbCost = getobjectivevalue(mp);
  ubCost = Inf;
  spbest = Dict();
  sqbest = Dict();

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
  sphatbar = copy(sphat);
  sqhatbar = copy(sqhat);
  # obtain the objective value to test whether the current solution is feasible
  vioSP = getobjectivevalue(sp);
  nonstop = true;
  while (nonstop)&(nC <= cMax)
    while (vioSP > ϵ*totalD)&(nC <= cMax)
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
      push!(C,cutTemp);

      # solve the master program, repeat the iteration
      # add the regularization term
      mpCopy = copy(mp);
      mpobj = mpCopy.obj;
      for i in fData.genIDList
        mpobj += ρ*((mpCopy[:sp][i] - sphatbar[i])^2 + (mpCopy[:sq][i] - sqhatbar[i])^2);
      end
      @objective(mpCopy,Min,mpobj);
      solve(mpCopy);
      # record the value of decision variable sp and sq
      sphat = Dict();
      sqhat = Dict();
      spobj = copy(spobjo);
      for i in fData.genIDList
        sphat[i] = getvalue(mpCopy[:sp][i]);
        sqhat[i] = getvalue(mpCopy[:sq][i]);
        # update the subproblem objective function with the current master solutions
        spobj -= (sphat[i]*sp[:λpi][fData.Loc[i]] + sqhat[i]*sp[:λqi][fData.Loc[i]]);
      end
      @objective(sp,Max,spobj);
      ubCost = costCal(fData,sphat,sqhat);
      solve(sp);
      vioSP = getobjectivevalue(sp);
      for i in fData.genIDList
        sphatbar[i] = sphat[i];
        sqhatbar[i] = sqhat[i];
      end
    end

    # compute the lower bound
    if ρ > 0
      solve(mp);
      if abs(getobjectivevalue(mp) - ubCost) < 1e-4*ubCost
        nonstop = false;
        lbCost = getobjectivevalue(mp);
        for i in fData.genIDList
          spbest[i] = sphatbar[i];
          sqbest[i] = sqhatbar[i];
        end
      else
        nonstop = true;
        spobj = copy(spobjo);
        for i in fData.genIDList
          sphat[i] = getvalue(mp[:sp][i]);
          sqhat[i] = getvalue(mp[:sq][i]);
          sphatbar[i] = sphat[i];
          sqhatbar[i] = sqhat[i];
          spobj -= (sphat[i]*sp[:λpi][fData.Loc[i]] + sqhat[i]*sp[:λqi][fData.Loc[i]]);
        end
        @objective(sp,Max,spobj);
        solve(sp);
        vioSP = getobjectivevalue(sp);
      end
    else
      nonstop = false;
      lbCost = getobjectivevalue(mp);
      for i in fData.genIDList
        spbest[i] = sphatbar[i];
        sqbest[i] = sqhatbar[i];
      end
      ubCost = costCal(fData,spbest,sqbest);
    end
  end

  solve(mp);
  lbCost = getobjectivevalue(mp);
  # sphat/sqhat is the final solution
  return spbest,sqbest,nC,C,vioSP,lbCost,ubCost;
end

function solveProcess_BudgetR_Corr(fData::fixedData,uData::Dict{Any,Any},vmax,vmin,θDmax,θDmin,ρ,partCorr,Γ = 1, ϵ = 1e-4, cMax = 500)
  # create the initial master/sub problem
  # without any scenario and cuts yet
  totalD = sum(abs(fData.Pd[j]) + abs(fData.Qd[j]) for j in fData.IDList);
  mp = createM(fData,uData);
  sp,spobjo = createS_Budget_CorrB(fData,uData,vmax,vmin,θDmax,θDmin,partCorr,Γ,numericOpt);
  # number of cuts so far and the set of cut coefficients
  nC = 0;
  C = [];

  solve(mp);
  lbCost = getobjectivevalue(mp);
  ubCost = Inf;
  spbest = Dict();
  sqbest = Dict();

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
  sphatbar = copy(sphat);
  sqhatbar = copy(sqhat);
  # obtain the objective value to test whether the current solution is feasible
  vioSP = getobjectivevalue(sp);
  nonstop = true;
  while (nonstop)&(nC <= cMax)
    while (vioSP > ϵ*totalD)&(nC <= cMax)
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
      push!(C,cutTemp);

      # solve the master program, repeat the iteration
      # add the regularization term
      mpCopy = copy(mp);
      mpobj = mpCopy.obj;
      for i in fData.genIDList
        mpobj += ρ*((mpCopy[:sp][i] - sphatbar[i])^2 + (mpCopy[:sq][i] - sqhatbar[i])^2);
      end
      @objective(mpCopy,Min,mpobj);
      solve(mpCopy);
      # record the value of decision variable sp and sq
      sphat = Dict();
      sqhat = Dict();
      spobj = copy(spobjo);
      for i in fData.genIDList
        sphat[i] = getvalue(mpCopy[:sp][i]);
        sqhat[i] = getvalue(mpCopy[:sq][i]);
        # update the subproblem objective function with the current master solutions
        spobj -= (sphat[i]*sp[:λpi][fData.Loc[i]] + sqhat[i]*sp[:λqi][fData.Loc[i]]);
      end
      @objective(sp,Max,spobj);
      ubCost = costCal(fData,sphat,sqhat);
      solve(sp);
      vioSP = getobjectivevalue(sp);
      for i in fData.genIDList
        sphatbar[i] = sphat[i];
        sqhatbar[i] = sqhat[i];
      end
    end

    # compute the lower bound
    if ρ > 0
      solve(mp);
      if abs(getobjectivevalue(mp) - ubCost) < 1e-4*ubCost
        nonstop = false;
        lbCost = getobjectivevalue(mp);
        for i in fData.genIDList
          spbest[i] = sphatbar[i];
          sqbest[i] = sqhatbar[i];
        end
      else
        nonstop = true;
        spobj = copy(spobjo);
        for i in fData.genIDList
          sphat[i] = getvalue(mp[:sp][i]);
          sqhat[i] = getvalue(mp[:sq][i]);
          sphatbar[i] = sphat[i];
          sqhatbar[i] = sqhat[i];
          spobj -= (sphat[i]*sp[:λpi][fData.Loc[i]] + sqhat[i]*sp[:λqi][fData.Loc[i]]);
        end
        @objective(sp,Max,spobj);
        solve(sp);
        vioSP = getobjectivevalue(sp);
      end
    else
      nonstop = false;
      lbCost = getobjectivevalue(mp);
      for i in fData.genIDList
        spbest[i] = sphatbar[i];
        sqbest[i] = sqhatbar[i];
      end
      ubCost = costCal(fData,spbest,sqbest);
    end
  end

  solve(mp);
  lbCost = getobjectivevalue(mp);
  # sphat/sqhat is the final solution
  return spbest,sqbest,nC,C,vioSP,lbCost,ubCost;
end
