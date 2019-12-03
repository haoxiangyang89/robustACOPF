# the extensive formulation of the rACOPF problem

function solveExtensiveNConv(fData,uData,Ω,yp,yq,vmax,vmin,θDmax,θDmin)
  # set up the extensive formulation
  mp = Model(solver = IpoptSolver(print_level = 0,max_iter = 10000,linear_solver = "ma27"));
  # obtain the pairs that are connected
  connectPair = [];
  connectDict = Dict();
  branchDict1 = Dict();
  branchDict2 = Dict();
  for i in fData.IDList
    connectDict[i] = [];
    branchDict1[i] = [];
    branchDict2[i] = [];
  end
  for k in fData.brList
    push!(branchDict1[k[1]],k);
    push!(branchDict2[k[2]],k);
    if !((k[1],k[2]) in connectPair)
      push!(connectPair,(k[1],k[2]));
      push!(connectDict[k[1]],k[2]);
    end
  end
  kpDict = Dict();
  for k in fData.brList
    if (k[1],k[2]) in keys(kpDict)
      push!(kpDict[(k[1],k[2])],k);
    else
      kpDict[(k[1],k[2])] = [k];
    end
  end

  con1Par = Dict();
  con2Par = Dict();
  for k in fData.brList
    if abs(fData.b[k]) > 10000
      con1Par[k] = 0.01;
      con2Par[k] = 100;
    else
      con1Par[k] = 1;
      con2Par[k] = 1;
    end
  end

  # set up the variables: active/reactive power injection
  @variable(mp, fData.Pmin[i] <= sp[i in fData.genIDList] <= fData.Pmax[i]);
  @variable(mp, fData.Qmin[i] <= sq[i in fData.genIDList] <= fData.Qmax[i]);

  sphatsum = Dict();
  for i in fData.IDList
    sphatsum[i] = @expression(mp,0.0);
    if i in keys(fData.LocRev)
      for j in fData.LocRev[i]
        sphatsum[i] += sp[j];
      end
    end
  end

  sqhatsum = Dict();
  for i in fData.IDList
    sqhatsum[i] = @expression(mp,0.0);
    if i in keys(fData.LocRev)
      for j in fData.LocRev[i]
        sqhatsum[i] += sq[j];
      end
    end
  end

  # obtain the bounds on uncertain injection for each bus
  upmax = Dict();
  upmin = Dict();
  uqmax = Dict();
  uqmin = Dict();
  for i in keys(uData)
    upmax[i] = uData[i].GPmax;
    upmin[i] = uData[i].GPmin;
    uqmax[i] = uData[i].GQmax;
    uqmin[i] = uData[i].GQmin;
  end

  LocRev = Dict();
  for i in fData.IDList
    if i in keys(fData.LocRev)
      LocRev[i] = fData.LocRev[i];
    else
      LocRev[i] = [];
    end
  end

  # set up the constraints:
  # upper bound/lower bound of the generator injection
  @constraint(mp, ubp[i in fData.IDList], sum(sp[j] for j in LocRev[i]) <= uData[i].rPmax - upmax[i]);
  @constraint(mp, ubq[i in fData.IDList], sum(sq[j] for j in LocRev[i]) <= uData[i].rQmax - uqmax[i]);
  @constraint(mp, lbp[i in fData.IDList], sum(sp[j] for j in LocRev[i]) >= uData[i].rPmin - upmin[i]);
  @constraint(mp, lbq[i in fData.IDList], sum(sq[j] for j in LocRev[i]) >= uData[i].rPmin - uqmin[i]);

  # set up the objective function
  @expression(mp,objExpr,0);
  for i in fData.genIDList
    cpGen = fData.cp[i];
    for j in 1:cpGen.n
      objExpr += cpGen.params[j]*(sp[i]*fData.baseMVA)^(cpGen.n-j);
    end
    if !(fData.cq == Dict())
      cqGen = fData.cq[i];
      for j in 1:cqGen.n
        objExpr += cqGen.params[j]*(sq[i]*fData.baseMVA)^(cqGen.n-j);
      end
    end
  end
  @objective(mp,Min,objExpr);

  @variable(mp,opplus[i in fData.IDList, ω in Ω] >= 0);
  @variable(mp,opminus[i in fData.IDList, ω in Ω] >= 0);
  @variable(mp,oqplus[i in fData.IDList, ω in Ω] >= 0);
  @variable(mp,oqminus[i in fData.IDList, ω in Ω] >= 0);

  @variable(mp,p[k in fData.brList, ω in Ω]);
  @variable(mp,ps[k in fData.brList, ω in Ω]);
  @variable(mp,q[k in fData.brList, ω in Ω]);
  @variable(mp,qs[k in fData.brList, ω in Ω]);
  @variable(mp,vmin[i] <= v[i in fData.IDList, ω in Ω] <= vmax[i]);
  @variable(mp,θ[i in fData.IDList, ω in Ω]);
  # set up the reference bus angle = 0
  for i in fData.IDList
    if fData.bType[i] == 3
      for ω in Ω
        @constraint(mp,θ[i,ω] == 0);
      end
    end
  end
  @constraint(mp,θDiff1[k in fData.brList, ω in Ω], (θ[k[1],ω] - fData.σ[k]) - θ[k[2],ω] <= θDmax[k]);
  @constraint(mp,θDiff2[k in fData.brList, ω in Ω], (θ[k[1],ω] - fData.σ[k]) - θ[k[2],ω] >= θDmin[k]);

  @NLconstraint(mp,pfp[k in fData.brList, ω in Ω], ps[k,ω] == con1Par[k]*(fData.g[k]*v[k[1],ω]^2/(fData.τ1[k]^2)
              - fData.g[k]*v[k[1],ω]*v[k[2],ω]*cos((θ[k[1],ω] - fData.σ[k]) - θ[k[2],ω])/(fData.τ1[k]*fData.τ2[k])
              - fData.b[k]*v[k[1],ω]*v[k[2],ω]*sin((θ[k[1],ω] - fData.σ[k]) - θ[k[2],ω])/(fData.τ1[k]*fData.τ2[k])));
  @NLconstraint(mp,qfp[k in fData.brList, ω in Ω], qs[k,ω] == con1Par[k]*((-fData.b[k] - fData.bc[k]/2)*v[k[1],ω]^2/(fData.τ1[k]^2)
              + fData.b[k]*v[k[1],ω]*v[k[2],ω]*cos((θ[k[1],ω] - fData.σ[k]) - θ[k[2],ω])/(fData.τ1[k]*fData.τ2[k])
              - fData.g[k]*v[k[1],ω]*v[k[2],ω]*sin((θ[k[1],ω] - fData.σ[k]) - θ[k[2],ω])/(fData.τ1[k]*fData.τ2[k])));
  @constraint(mp,pScale[k in fData.brList, ω in Ω], p[k,ω] == con2Par[k]*ps[k,ω]);
  @constraint(mp,qScale[k in fData.brList, ω in Ω], q[k,ω] == con2Par[k]*qs[k,ω]);
  @NLconstraint(mp,flowConstr[k in fData.brList, ω in Ω], p[k,ω]^2 + q[k,ω]^2 <= fData.rateA[k]^2);

  # set up the flow balance
  @constraint(mp,pbalance[i in fData.IDList,ω in Ω],sum(p[k,ω] for k in branchDict1[i]) + v[i,ω]^2*fData.gs[i] + opplus[i,ω] - opminus[i,ω] ==
    sphatsum[i] + uData[i].GP0 + (uData[i].GPmax - uData[i].GP0)*(max(2*yp[ω,fData.busInd[i]] - 1,0)) + (uData[i].GPmin - uData[i].GP0)*(max(1 - 2*yp[ω,fData.busInd[i]],0)));
  @constraint(mp,qbalance[i in fData.IDList,ω in Ω],sum(q[k,ω] for k in branchDict1[i]) - v[i,ω]^2*fData.bs[i] + oqplus[i,ω] - oqminus[i,ω] ==
    sqhatsum[i] + uData[i].GQ0 + (uData[i].GQmax - uData[i].GQ0)*(max(2*yq[ω,fData.busInd[i]] - 1,0)) + (uData[i].GQmin - uData[i].GQ0)*(max(1 - 2*yq[ω,fData.busInd[i]],0)));

  @constraint(mp,opConstr1[i in fData.IDList, ω in Ω],opminus[i,ω] <= uData[i].WPmax);
  @constraint(mp,oqConstr1[i in fData.IDList, ω in Ω],oqminus[i,ω] <= uData[i].WQmax);
  @constraint(mp,opConstr2[i in fData.IDList, ω in Ω],opplus[i,ω] <= uData[i].WPmax + uData[i].RESP0 +
    (uData[i].RESPmax - uData[i].RESP0)*(max(2*yp[ω,fData.busInd[i]] - 1,0)) + (uData[i].RESPmin - uData[i].RESP0)*(max(1 - 2*yp[ω,fData.busInd[i]],0)));
  @constraint(mp,oqConstr2[i in fData.IDList, ω in Ω],oqplus[i,ω] <= uData[i].WQmax + uData[i].RESQ0 +
    (uData[i].RESQmax - uData[i].RESQ0)*(max(2*yq[ω,fData.busInd[i]] - 1,0)) + (uData[i].RESQmin - uData[i].RESQ0)*(max(1 - 2*yq[ω,fData.busInd[i]],0)));
  return mp;
end

# an iterative process to append important scenarios to solve the extensive nonconvex case
function solveENC_Proc(fData,uData,Ωs,yps,yqs,vmax,vmin,θDmax,θDmin,ϵ = 1e-4)
    # initialze with the first scenario
    totalD = sum(abs(fData.Pd[j]) + abs(fData.Qd[j]) for j in fData.IDList);
    fixedy = [1];
    contBool = true;
    sphatnc = Dict();
    sqhatnc = Dict();
    ncnbC = Inf;
    vio = Inf;

    # test if the solution is nonconvex robust feasible
    while contBool
        mext = solveExtensiveNConv(fData,uData,1:length(fixedy),yps[fixedy,:],yqs[fixedy,:],vmax,vmin,θDmax,θDmin);
        meStatus = solve(mext);
        if meStatus == :Infeasible
            ncnbC = Inf;
            contBool = false;
        else
            ncnbC = getobjectivevalue(mext);
            for i in fData.genIDList
                sphatnc[i] = getvalue(mext[:sp][i]);
                sqhatnc[i] = getvalue(mext[:sq][i]);
            end
            vioNCListb = testSolNC(fData,uData,sphatnc,sqhatnc,Ωs,yps,yqs,"L",vmax,vmin,θDmax,θDmin);
            if (maximum(vioNCListb) > ϵ*totalD)&&(length(fixedy) <= 20)
                push!(fixedy,indmax(vioNCListb));
                println(ncnbC," ",maximum(vioNCListb));
                if vio > maximum(vioNCListb)
                  vio = maximum(vioNCListb);
                end
            else
                contBool = false;
                vio = maximum(vioNCListb);
            end
        end
    end
    return sphatnc,sqhatnc,ncnbC,vio,length(fixedy);
end

# solve the extensive fomulation of the convexly relaxed rACOPF problem
function solveExtensiveConv(fData,uData,Ω,yp,yq,vmax,vmin,θDmax,θDmin,βSwitch = 1)
  θu = Dict();
  for k in fData.brList
    θu[k] = max(abs(θDmax[k]),abs(θDmin[k]));
  end
  # scaling issue exists for this primal problem.
  mp = Model(solver = IpoptSolver(print_level = 0,max_iter = 10000,linear_solver = "ma27"));

  # obtain the pairs that are connected
  connectPair = [];
  connectDict = Dict();
  branchDict1 = Dict();
  branchDict2 = Dict();
  for i in fData.IDList
    connectDict[i] = [];
    branchDict1[i] = [];
    branchDict2[i] = [];
  end
  for k in fData.brList
    push!(branchDict1[k[1]],k);
    push!(branchDict2[k[2]],k);
    if !((k[1],k[2]) in connectPair)
      push!(connectPair,(k[1],k[2]));
      push!(connectDict[k[1]],k[2]);
    end
  end
  kpDict = Dict();
  for k in fData.brList
    if (k[1],k[2]) in keys(kpDict)
      push!(kpDict[(k[1],k[2])],k);
    else
      kpDict[(k[1],k[2])] = [k];
    end
  end

  # set up the variables: active/reactive power injection
  @variable(mp, fData.Pmin[i] <= sp[i in fData.genIDList] <= fData.Pmax[i]);
  @variable(mp, fData.Qmin[i] <= sq[i in fData.genIDList] <= fData.Qmax[i]);
  sphatsum = Dict();
  for i in fData.IDList
    sphatsum[i] = @expression(mp,0.0);
    if i in keys(fData.LocRev)
      for j in fData.LocRev[i]
        sphatsum[i] += sp[j];
      end
    end
  end

  sqhatsum = Dict();
  for i in fData.IDList
    sqhatsum[i] = @expression(mp,0.0);
    if i in keys(fData.LocRev)
      for j in fData.LocRev[i]
        sqhatsum[i] += sq[j];
      end
    end
  end

  # set up the objective function
  @expression(mp,objExpr,0);
  for i in fData.genIDList
    cpGen = fData.cp[i];
    for j in 1:cpGen.n
      objExpr += cpGen.params[j]*(sp[i]*fData.baseMVA)^(cpGen.n-j);
    end
    if !(fData.cq == Dict())
      cqGen = fData.cq[i];
      for j in 1:cqGen.n
        objExpr += cqGen.params[j]*(sq[i]*fData.baseMVA)^(cqGen.n-j);
      end
    end
  end
  @objective(mp,Min,objExpr);

  @variable(mp,p[k in fData.brList, ω in Ω]);
  @variable(mp,q[k in fData.brList, ω in Ω]);
  @variable(mp,ps[k in fData.brList, ω in Ω]);
  @variable(mp,qs[k in fData.brList, ω in Ω]);
  @variable(mp,vmin[i] <= v[i in fData.IDList, ω in Ω] <= vmax[i]);
  @variable(mp,vhat[i in fData.IDList, ω in Ω]);
  @variable(mp,θ[i in fData.IDList, ω in Ω]);
  # set up the reference bus angle = 0
  for i in fData.IDList
    if fData.bType[i] == 3
      for ω in Ω
        @constraint(mp,θ[i,ω] == 0);
      end
    end
  end

  @variable(mp,opplus[i in fData.IDList, ω in Ω] >= 0);
  @variable(mp,opminus[i in fData.IDList, ω in Ω] >= 0);
  @variable(mp,oqplus[i in fData.IDList, ω in Ω] >= 0);
  @variable(mp,oqminus[i in fData.IDList, ω in Ω] >= 0);

  @variable(mp,vv[k in fData.brList, ω in Ω]);
  @variable(mp,cos(θu[k]) <= cs[k in fData.brList, ω in Ω] <= 1);
  @variable(mp,-sin(θu[k]) <= ss[k in fData.brList, ω in Ω] <= sin(θu[k]));
  # add the strengthened bounds on cs
  csmax = Dict();
  csmin = Dict();
  ssmax = Dict();
  ssmin = Dict();
  for k in fData.brList
    for ω in Ω
      if (θDmax[k] >= 0)&(θDmin[k] >= 0)
        csmax[k] = cos(θDmin[k]);
        csmin[k] = cos(θDmax[k]);
        @constraint(mp, cs[k,ω] <= csmax[k]);
        @constraint(mp, cs[k,ω] >= csmin[k]);
      elseif (θDmax[k] < 0)&(θDmin[k] < 0)
        csmax[k] = cos(θDmax[k]);
        csmin[k] = cos(θDmin[k]);
        @constraint(mp, cs[k,ω] <= csmax[k]);
        @constraint(mp, cs[k,ω] >= csmin[k]);
      else
        csmax[k] = 1;
        csmin[k] = min(cos(θDmax[k]),cos(θDmin[k]));
        @constraint(mp, cs[k,ω] <= csmax[k]);
        @constraint(mp, cs[k,ω] >= csmin[k]);
      end
      @constraint(mp, cs[k,ω] >= (cos(θDmax[k]) - cos(θDmin[k]))/(θDmax[k] - θDmin[k])*((θ[k[1],ω] - fData.σ[k]) - θ[k[2],ω] - θDmin[k]) + cos(θDmin[k]));
    end
  end
  # add the strengthened bounds on ss
  for k in fData.brList
    for ω in Ω
      ssmax[k] = sin(θDmax[k]);
      ssmin[k] = sin(θDmin[k]);
      @constraint(mp, ss[k,ω] <= ssmax[k]);
      @constraint(mp, ss[k,ω] >= ssmin[k]);
    end
  end

  @variable(mp,wc[k in fData.brList, ω in Ω]);
  @variable(mp,ws[k in fData.brList, ω in Ω]);
  @constraint(mp,wcEquality[k in fData.brList, ω in Ω;k[1] < k[2]], wc[k,ω] == wc[(k[2],k[1],k[3]),ω]);
  @constraint(mp,wsEquality[k in fData.brList, ω in Ω;k[1] < k[2]], ws[k,ω] == -ws[(k[2],k[1],k[3]),ω]);
  @constraint(mp,vvEquality[k in fData.brList, ω in Ω;k[1] < k[2]], vv[k,ω] == vv[(k[2],k[1],k[3]),ω]);
  @constraint(mp,csEquality[k in fData.brList, ω in Ω;k[1] < k[2]], cs[k,ω] == cs[(k[2],k[1],k[3]),ω]);
  @constraint(mp,ssEquality[k in fData.brList, ω in Ω;k[1] < k[2]], ss[k,ω] == -ss[(k[2],k[1],k[3]),ω]);

  con1Par = Dict();
  con2Par = Dict();
  for k in fData.brList
    con1Par[k] = 1/sqrt(abs(fData.b[k]));
    con2Par[k] = sqrt(abs(fData.b[k]));
  end

  @constraint(mp,lineConstrP[k in fData.brList, ω in Ω],ps[k,ω] == con1Par[k]*(fData.g[k]*vhat[k[1],ω]/(fData.τ1[k]^2) - fData.g[k]*wc[k,ω] - fData.b[k]*ws[k,ω]));
  @constraint(mp,pScale[k in fData.brList, ω in Ω],p[k,ω] == con2Par[k]*ps[k,ω]);
  @constraint(mp,lineConstrQ[k in fData.brList, ω in Ω],qs[k,ω] == con1Par[k]*((-fData.b[k] - fData.bc[k]/2)*vhat[k[1],ω]/(fData.τ1[k]^2) + fData.b[k]*wc[k,ω] - fData.g[k]*ws[k,ω]));
  @constraint(mp,qScale[k in fData.brList, ω in Ω],q[k,ω] == con2Par[k]*qs[k,ω]);

  @constraint(mp,socConstraint1[k in fData.brList, ω in Ω; fData.rateA[k] < Inf],p[k,ω]^2 + q[k,ω]^2 <= fData.rateA[k]^2);
  @NLconstraint(mp,socConstraint2[k in fData.brList, ω in Ω],vv[k,ω]^2 <= vhat[k[1],ω]/fData.τ1[k]^2 * vhat[k[2],ω]/fData.τ2[k]^2);
  @NLconstraint(mp,socConstraint3[k in fData.brList, ω in Ω],cs[k,ω] + (1-cos(θu[k]))/(θu[k])^2*((θ[k[1],ω] - fData.σ[k]) - θ[k[2],ω])^2 <= 1);
  @NLconstraint(mp,socConstraint4[i in fData.IDList, ω in Ω],v[i,ω]^2 <= vhat[i,ω]);
  @NLconstraint(mp,socConstraint5[k in fData.brList, ω in Ω],wc[k,ω]^2 + ws[k,ω]^2 <= vhat[k[1],ω]/(fData.τ1[k]^2)*vhat[k[2],ω]/(fData.τ2[k]^2));

  @constraint(mp,sinMock1[k in fData.brList, ω in Ω],
              ss[k,ω] <= cos(θu[k]/2)*((θ[k[1],ω] - fData.σ[k]) - θ[k[2],ω] - θu[k]/2) + sin(θu[k]/2));
  @constraint(mp,sinMock2[k in fData.brList, ω in Ω],
              ss[k,ω] >= cos(θu[k]/2)*((θ[k[1],ω] - fData.σ[k]) - θ[k[2],ω] + θu[k]/2) - sin(θu[k]/2));
  @constraint(mp,angleDiff1[k in fData.brList, ω in Ω],(θ[k[1],ω] - fData.σ[k]) - θ[k[2],ω] <= θDmax[k]);
  @constraint(mp,angleDiff2[k in fData.brList, ω in Ω],(θ[k[1],ω] - fData.σ[k]) - θ[k[2],ω] >= θDmin[k]);
  @constraint(mp,v2Mock2[i in fData.IDList, ω in Ω], vhat[i,ω] - (vmax[i] + vmin[i])*v[i,ω] <= -vmax[i]*vmin[i]);
  @constraint(mp,vvMock1[k in fData.brList, ω in Ω],
              vv[k,ω] >= vmin[k[1]]*v[k[2],ω]/(fData.τ1[k]*fData.τ2[k]) + vmin[k[2]]*v[k[1],ω]/(fData.τ1[k]*fData.τ2[k]) - vmin[k[1]]*vmin[k[2]]/(fData.τ1[k]*fData.τ2[k]));
  @constraint(mp,vvMock2[k in fData.brList, ω in Ω],
              vv[k,ω] >= vmax[k[1]]*v[k[2],ω]/(fData.τ1[k]*fData.τ2[k]) + vmax[k[2]]*v[k[1],ω]/(fData.τ1[k]*fData.τ2[k]) - vmax[k[1]]*vmax[k[2]]/(fData.τ1[k]*fData.τ2[k]));
  @constraint(mp,vvMock3[k in fData.brList, ω in Ω],
              vv[k,ω] <= vmin[k[1]]*v[k[2],ω]/(fData.τ1[k]*fData.τ2[k]) + vmax[k[2]]*v[k[1],ω]/(fData.τ1[k]*fData.τ2[k]) - vmin[k[1]]*vmax[k[2]]/(fData.τ1[k]*fData.τ2[k]));
  @constraint(mp,vvMock4[k in fData.brList, ω in Ω],
              vv[k,ω] <= vmax[k[1]]*v[k[2],ω]/(fData.τ1[k]*fData.τ2[k]) + vmin[k[2]]*v[k[1],ω]/(fData.τ1[k]*fData.τ2[k]) - vmax[k[1]]*vmin[k[2]]/(fData.τ1[k]*fData.τ2[k]));

  @constraint(mp,wcMock1[k in fData.brList, ω in Ω],
              wc[k,ω] >= vmin[k[1]]*vmin[k[2]]/(fData.τ1[k]*fData.τ2[k])*cs[k,ω] + csmin[k]*vv[k,ω] - vmin[k[1]]*vmin[k[2]]/(fData.τ1[k]*fData.τ2[k])*csmin[k]);
  @constraint(mp,wcMock2[k in fData.brList, ω in Ω],
              wc[k,ω] >= vmax[k[1]]*vmax[k[2]]/(fData.τ1[k]*fData.τ2[k])*cs[k,ω] + csmax[k]*vv[k,ω] - vmax[k[1]]*vmax[k[2]]/(fData.τ1[k]*fData.τ2[k])*csmax[k]);
  @constraint(mp,wcMock3[k in fData.brList, ω in Ω],
              wc[k,ω] <= vmin[k[1]]*vmin[k[2]]/(fData.τ1[k]*fData.τ2[k])*cs[k,ω] + vv[k,ω]*csmax[k] - vmin[k[1]]*vmin[k[2]]/(fData.τ1[k]*fData.τ2[k])*csmax[k]);
  @constraint(mp,wcMock4[k in fData.brList, ω in Ω],
              wc[k,ω] <= vmax[k[1]]*vmax[k[2]]/(fData.τ1[k]*fData.τ2[k])*cs[k,ω] + vv[k,ω]*csmin[k] - vmax[k[1]]*vmax[k[2]]/(fData.τ1[k]*fData.τ2[k])*csmin[k]);

  @constraint(mp,wsMock1[k in fData.brList, ω in Ω],
              ws[k,ω] >= vmin[k[1]]*vmin[k[2]]/(fData.τ1[k]*fData.τ2[k])*ss[k,ω] + ssmin[k]*vv[k,ω] - vmin[k[1]]*vmin[k[2]]/(fData.τ1[k]*fData.τ2[k])*ssmin[k]);
  @constraint(mp,wsMock2[k in fData.brList, ω in Ω],
              ws[k,ω] >= vmax[k[1]]*vmax[k[2]]/(fData.τ1[k]*fData.τ2[k])*ss[k,ω] + ssmax[k]*vv[k,ω] - vmax[k[1]]*vmax[k[2]]/(fData.τ1[k]*fData.τ2[k])*ssmax[k]);
  @constraint(mp,wsMock3[k in fData.brList, ω in Ω],
              ws[k,ω] <= vmin[k[1]]*vmin[k[2]]/(fData.τ1[k]*fData.τ2[k])*ss[k,ω] + vv[k,ω]*ssmax[k] - vmin[k[1]]*vmin[k[2]]/(fData.τ1[k]*fData.τ2[k])*ssmax[k]);
  @constraint(mp,wsMock4[k in fData.brList, ω in Ω],
              ws[k,ω] <= vmax[k[1]]*vmax[k[2]]/(fData.τ1[k]*fData.τ2[k])*ss[k,ω] + vv[k,ω]*ssmin[k] - vmax[k[1]]*vmax[k[2]]/(fData.τ1[k]*fData.τ2[k])*ssmin[k]);

  @constraint(mp,tanConstr1[k in fData.brList, ω in Ω],
                ws[k,ω] - tan(θDmax[k])*wc[k,ω] <= 0);
  @constraint(mp,tanConstr2[k in fData.brList, ω in Ω],
                ws[k,ω] - tan(θDmin[k])*wc[k,ω] >= 0);

  vδ = Dict();
  θϕ = Dict();
  θδ = Dict();
  for i in fData.IDList
    vδ[i] = vmax[i]+vmin[i];
  end
  for k in fData.brList
    θϕ[k] = (θDmax[k] + θDmin[k])/2;
    θδ[k] = (θDmax[k] - θDmin[k])/2;
  end

  @constraint(mp,lncConstr1[k in fData.brList, ω in Ω],
              vδ[k[1]]*vδ[k[2]]/(fData.τ1[k]*fData.τ2[k])*(wc[k,ω]*cos(θϕ[k]) + ws[k,ω]*sin(θϕ[k])) -
              vmax[k[2]]/fData.τ2[k]*cos(θδ[k])*vδ[k[2]]/fData.τ2[k]*vhat[k[1],ω]/(fData.τ1[k]^2) -
              vmax[k[1]]/fData.τ1[k]*cos(θδ[k])*vδ[k[1]]/fData.τ1[k]*vhat[k[2],ω]/(fData.τ2[k]^2) >=
              vmax[k[1]]*vmax[k[2]]/(fData.τ1[k]*fData.τ2[k])*cos(θδ[k])*(vmin[k[1]]*vmin[k[2]]/(fData.τ1[k]*fData.τ2[k]) -
              vmax[k[1]]*vmax[k[2]]/(fData.τ1[k]*fData.τ2[k])));
  @constraint(mp,lncConstr2[k in fData.brList, ω in Ω],
              vδ[k[1]]*vδ[k[2]]/(fData.τ1[k]*fData.τ2[k])*(wc[k,ω]*cos(θϕ[k]) + ws[k,ω]*sin(θϕ[k])) -
              vmin[k[2]]/fData.τ2[k]*cos(θδ[k])*vδ[k[2]]/fData.τ2[k]*vhat[k[1],ω]/(fData.τ1[k]^2) -
              vmin[k[1]]/fData.τ1[k]*cos(θδ[k])*vδ[k[1]]/fData.τ1[k]*vhat[k[2],ω]/(fData.τ2[k]^2) >=
              -vmin[k[1]]*vmin[k[2]]/(fData.τ1[k]*fData.τ2[k])*cos(θδ[k])*(vmin[k[1]]*vmin[k[2]]/(fData.τ1[k]*fData.τ2[k]) -
              vmax[k[1]]*vmax[k[2]]/(fData.τ1[k]*fData.τ2[k])));

  @constraint(mp,opConstr1[i in fData.IDList, ω in Ω],opminus[i,ω] <= uData[i].WPmax);
  @constraint(mp,oqConstr1[i in fData.IDList, ω in Ω],oqminus[i,ω] <= uData[i].WQmax);
  @constraint(mp,opConstr2[i in fData.IDList, ω in Ω],opplus[i,ω] <= uData[i].WPmax + uData[i].RESP0 +
    (uData[i].RESPmax - uData[i].RESP0)*(max(2*yp[ω,fData.busInd[i]] - 1,0)) + (uData[i].RESPmin - uData[i].RESP0)*(max(1 - 2*yp[ω,fData.busInd[i]],0)));
  @constraint(mp,oqConstr2[i in fData.IDList, ω in Ω],oqplus[i,ω] <= uData[i].WQmax + uData[i].RESQ0 +
    (uData[i].RESQmax - uData[i].RESQ0)*(max(2*yq[ω,fData.busInd[i]] - 1,0)) + (uData[i].RESQmin - uData[i].RESQ0)*(max(1 - 2*yq[ω,fData.busInd[i]],0)));

  @constraint(mp,pbalance[i in fData.IDList,ω in Ω],sum(p[k,ω] for k in branchDict1[i]) + vhat[i,ω]*fData.gs[i] + opplus[i,ω] - opminus[i,ω] ==
    sphatsum[i] + uData[i].GP0 + (uData[i].GPmax - uData[i].GP0)*(max(2*yp[ω,fData.busInd[i]] - 1,0)) + (uData[i].GPmin - uData[i].GP0)*(max(1 - 2*yp[ω,fData.busInd[i]],0)));
  @constraint(mp,qbalance[i in fData.IDList,ω in Ω],sum(q[k,ω] for k in branchDict1[i]) - vhat[i,ω]*fData.bs[i] + oqplus[i,ω] - oqminus[i,ω] ==
    sqhatsum[i] + uData[i].GQ0 + (uData[i].GQmax - uData[i].GQ0)*(max(2*yq[ω,fData.busInd[i]] - 1,0)) + (uData[i].GQmin - uData[i].GQ0)*(max(1 - 2*yq[ω,fData.busInd[i]],0)));

  return mp;
end

function solveEC_Proc(fData,uData,Ωs,yps,yqs,vmax,vmin,θDmax,θDmin,ϵ = 1e-4)
    # initialze with the first scenario
    totalD = sum(abs(fData.Pd[j]) + abs(fData.Qd[j]) for j in fData.IDList);
    fixedy = [1];
    contBool = true;
    sphatb = Dict();
    sqhatb = Dict();
    cbC = Inf;
    vio = Inf;

    # test if the solution is nonconvex robust feasible
    while contBool
        mext = solveExtensiveConv(fData,uData,1:length(fixedy),yps[fixedy,:],yqs[fixedy,:],vmax,vmin,θDmax,θDmin);
        meStatus = solve(mext);
        if meStatus == :Infeasible
            cbC = Inf;
            contBool = false;
        else
            cbC = getobjectivevalue(mext);
            for i in fData.genIDList
                sphatb[i] = getvalue(mext[:sp][i]);
                sqhatb[i] = getvalue(mext[:sq][i]);
            end
            vioCListb = testSolC(fData,uData,sphatb,sqhatb,Ωs,yps,yqs,"L",vmax,vmin,θDmax,θDmin);
            if (maximum(vioCListb) > ϵ*totalD)&&(length(fixedy) <= 20)
                push!(fixedy,indmax(vioCListb));
                println(cbC," ",maximum(vioCListb));
                if vio > maximum(vioCListb)
                  vio = maximum(vioCListb);
                end
            else
                contBool = false;
                vio = maximum(vioCListb);
            end
        end
    end
    return sphatb,sqhatb,cbC,vio,length(fixedy);
end
