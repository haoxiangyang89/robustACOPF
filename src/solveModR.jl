# this is the function to execute the solving process
# the output will be the optimal solution
function solveProcessR(fData::fixedData,uData::Dict{Any,Any},vmax,vmin,θDmax,θDmin,ρ,ϵ = 1e-4)
  # create the initial master/sub problem
  # without any scenario and cuts yet
  mp = createM(fData,uData);
  sp,spobjo = createS(fData,uData,vmax,vmin,θDmax,θDmin);
  # number of cuts so far and the set of cut coefficients
  nC = 0;
  C = [];

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
#  sp.solver = GurobiSolver(MIPGap = 0.04,MIQCPMethod = 0);
  solve(sp);
  # obtain the objective value to test whether the current solution is feasible
  vioSP = getobjectivevalue(sp);
  while vioSP > ϵ
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
    spobj = copy(spobjo);
    for i in fData.genIDList
      sphat[i] = getvalue(mpCopy[:sp][i]);
      sqhat[i] = getvalue(mpCopy[:sq][i]);
      # update the subproblem objective function with the current master solutions
      spobj -= (sphat[i]*sp[:λpi][fData.Loc[i]] + sqhat[i]*sp[:λqi][fData.Loc[i]]);
    end
    @objective(sp,Max,spobj);
    solve(sp);
    # obtain the objective value to test whether the current solution is feasible
    vioSP = getobjectivevalue(sp);
  end

  ypr = Dict();
  yqr = Dict();
  for i in fData.IDList
    ypr[i] = getvalue(sp[:yp][i]);
  end
  for i in fData.IDList
    yqr[i] = getvalue(sp[:yq][i]);
  end

  solve(mp);
  lbCost = getobjectivevalue(mp);
  # sphat/sqhat is the final solution
  return sphat,sqhat,nC,C,vioSP,ypr,yqr,lbCost;
end
