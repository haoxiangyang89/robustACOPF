# Nonconvex ACOPF Nominal Formulation
function ACOPF_NC(fData,uData,vmax,vmin,θDmax,θDmin,βSwitch = 1)
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

  mp = Model(solver = IpoptSolver(print_level = 0,max_iter = 10000,linear_solver = "ma27"));
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

  LocRev = Dict();
  for i in fData.IDList
    if i in keys(fData.LocRev)
      LocRev[i] = fData.LocRev[i];
    else
      LocRev[i] = [];
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

  @variable(mp,p[k in fData.brList]);
  @variable(mp,q[k in fData.brList]);
  @variable(mp,vmin[i] <= v[i in fData.IDList] <= vmax[i]);
  @variable(mp,θ[i in fData.IDList]);

  @variable(mp,opplus[i in fData.IDList] >= 0);
  @variable(mp,opminus[i in fData.IDList] >= 0);
  @variable(mp,oqplus[i in fData.IDList] >= 0);
  @variable(mp,oqminus[i in fData.IDList] >= 0);

  @constraint(mp,θDiff1[k in fData.brList], θ[k[1]] - fData.σ[k] - θ[k[2]] <= θDmax[k]);
  @constraint(mp,θDiff2[k in fData.brList], θ[k[1]] - fData.σ[k] - θ[k[2]] >= θDmin[k]);

  @NLconstraint(mp,pfp[k in fData.brList], p[k] == fData.g[k]*v[k[1]]^2/(fData.τ1[k]^2)
              - fData.g[k]*(v[k[1]]/fData.τ1[k])*(v[k[2]]/fData.τ2[k])*cos(θ[k[1]] - fData.σ[k] - θ[k[2]])
              - fData.b[k]*(v[k[1]]/fData.τ1[k])*(v[k[2]]/fData.τ2[k])*sin(θ[k[1]] - fData.σ[k] - θ[k[2]]));
  @NLconstraint(mp,qfp[k in fData.brList], q[k] == (-fData.b[k] - fData.bc[k]/2)*v[k[1]]^2/(fData.τ1[k]^2)
              + fData.b[k]*(v[k[1]]/fData.τ1[k])*(v[k[2]]/fData.τ2[k])*cos(θ[k[1]] - fData.σ[k] - θ[k[2]])
              - fData.g[k]*(v[k[1]]/fData.τ1[k])*(v[k[2]]/fData.τ2[k])*sin(θ[k[1]] - fData.σ[k] - θ[k[2]]));
  @NLconstraint(mp,flowConstr[k in fData.brList], p[k]^2 + q[k]^2 <= fData.rateA[k]^2);

  conWPpar1 = Dict();
  conWQpar1 = Dict();
  conWPpar2 = Dict();
  conWQpar2 = Dict();
  for i in fData.IDList
    if (abs(uData[i].WPmax) <= 1e-3)&(abs(uData[i].WPmax) < 0)
      conWPpar1[i] = 1/sqrt(abs(uData[i].WPmax));
      conWPpar2[i] = sqrt(abs(uData[i].WPmax));
    else
      conWPpar1[i] = 1;
      conWPpar2[i] = 1;
    end
    if (abs(uData[i].WQmax) <= 1e-3)&(abs(uData[i].WQmax) < 0)
      conWQpar1[i] = 1/sqrt(abs(uData[i].WQmax));
      conWQpar2[i] = sqrt(abs(uData[i].WQmax));
    else
      conWQpar1[i] = 1;
      conWQpar2[i] = 1;
    end
  end
  if βSwitch == 1
    @constraint(mp,opConstr1[i in fData.IDList],opminus[i] <= conWQpar1[i]*uData[i].WPmax);
    @constraint(mp,oqConstr1[i in fData.IDList],oqminus[i] <= conWQpar1[i]*uData[i].WQmax);
    @constraint(mp,opConstr2[i in fData.IDList],opplus[i] <= conWQpar1[i]*(uData[i].WPmax + uData[i].RESP0));
    @constraint(mp,oqConstr2[i in fData.IDList],oqplus[i] <= conWQpar1[i]*(uData[i].WQmax + uData[i].RESQ0));
  else
    @constraint(mp,opConstr[i in fData.IDList],opplus[i] + opminus[i] <= 0);
    @constraint(mp,oqConstr[i in fData.IDList],oqplus[i] + oqminus[i] <= 0);
  end

  @constraint(mp,pbalance[i in fData.IDList],sum(p[k] for k in branchDict1[i]) + v[i]^2*fData.gs[i]
      + opplus[i]*conWPpar2[i] - opminus[i]*conWPpar2[i] == sphatsum[i] + uData[i].GP0);
  @constraint(mp,qbalance[i in fData.IDList],sum(q[k] for k in branchDict1[i]) - v[i]^2*fData.bs[i]
      + oqplus[i]*conWQpar2[i] - oqminus[i]*conWQpar2[i] == sqhatsum[i] + uData[i].GQ0);

  return mp;
end

function ACOPF_C(fData,uData,vmax,vmin,θDmax,θDmin,βSwitch = 1)
  θu = Dict();
  for k in fData.brList
    θu[k] = max(abs(θDmax[k]),abs(θDmin[k]));
  end
  # scaling issue exists for this primal problem.
  mp = Model(solver = GurobiSolver(NumericFocus = 3, Threads = 20));
  # mp = Model(solver = MosekSolver(MSK_DPAR_MIO_TOL_REL_GAP = 1e-7,MSK_DPAR_INTPNT_TOL_PFEAS = 1e-9,MSK_DPAR_MIO_TOL_FEAS = 1e-9,MSK_IPAR_LOG = 0));

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
  @variable(mp, sps[i in fData.genIDList]);
  @variable(mp, fData.Qmin[i] <= sq[i in fData.genIDList] <= fData.Qmax[i]);
  @variable(mp, sqs[i in fData.genIDList]);
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

  # construct the objective function
  @expression(mp,objExpr,0);
  xp = Dict();
  xq = Dict();
  zp = Dict();
  zq = Dict();
  @variable(mp,tp[i in fData.genIDList]);
  @variable(mp,tq[i in fData.genIDList]);
  objPar1p = Dict();
  objPar2p = Dict();
  objPar1q = Dict();
  objPar2q = Dict();
  for i in fData.genIDList
    cpGen = fData.cp[i];
    if (cpGen.model == 1)
      if cpGen.n > 4
        # use anonymous constraints here since we don't really need to know the cost calculation procedure
        xp[i] = @variable(mp,[1:(div(cpGen.n,2)-1)],Bin);
        zp[i] = @variable(mp);
        @constraint(mp,sp[i] <= cpGen.params[cpGen.n-1]);
        @constraint(mp,sp[i] >= cpGen.params[1]);
        @constraint(mp,sum(xp[i][j] for j in 1:(div(cpGen.n,2)-1)) == 1);
        @constraint(mp,[j in 1:(div(cpGen.n,2)-1)],sp[i] <= (1-xp[i][j])*(cpGen.params[cpGen.n-1]-cpGen.params[1]+10)+cpGen.params[(j+1)*2]);
        @constraint(mp,[j in 1:(div(cpGen.n,2)-1)],sp[i] >= (xp[i][j]-1)*(cpGen.params[cpGen.n-1]-cpGen.params[1]+10)+cpGen.params[j*2]);
        @constraint(mp,[j in 1:(div(cpGen.n,2)-1)],zp[i] >= sp[i]*(cpGen.params[(j+1)*2] - cpGen.params[j*2])/
          (cpGen.params[(j+1)*2-1] - cpGen.params[j*2-1])+(cpGen.params[j*2]-(cpGen.params[(j+1)*2] - cpGen.params[j*2])/
          (cpGen.params[(j+1)*2-1] - cpGen.params[j*2-1])*cpGen.params[j*2-1]));
        objExpr += zp[j];
      elseif cpGen.n == 4
        @constraint(mp,sp[i] <= cpGen.params[cpGen.n-1]);
        @constraint(mp,sp[i] >= cpGen.params[1]);
        objExpr += (cpGen.params[4] - cpGen.params[2])/(cpGen.params[3] - cpGen.params[1])*sp[i] + cpGen.params[2];
      else
        error("Not enough parameters");
      end
    else
      if (cpGen.n == 3)&(cpGen.params[1] > 0)
        # quadratic
        αi = sqrt(cpGen.params[1])*fData.baseMVA;
        βi = cpGen.params[2]/(2*sqrt(cpGen.params[1]));
        f1i = -1/4 + cpGen.params[2]^2/(4*cpGen.params[1]) - cpGen.params[3];
        f2i = 1/4 + cpGen.params[2]^2/(4*cpGen.params[1]) - cpGen.params[3];
        @constraint(mp,norm([αi*sp[i] + βi, tp[i] + f1i]) <= tp[i] + f2i);
        @constraint(mp,tp[i] + f2i >= 0);
      elseif (cpGen.n == 3)&(cpGen.params[1] == 0)
        # quadratic but the quadratic coefficient is 0, so linear
        if cpGen.params[2] > 0
          objPar1p[i] = 1/sqrt(cpGen.params[2]*fData.baseMVA);
          objPar2p[i] = sqrt(cpGen.params[2]*fData.baseMVA);
        else
          objPar1p[i] = 1;
          objPar2p[i] = 1;
        end
        @constraint(mp,tp[i] == objPar1p[i]*(cpGen.params[2]*fData.baseMVA*sps[i]) + cpGen.params[3]);
        @constraint(mp,sps[i] == sp[i]*objPar2p[i]);
      elseif (cpGen.n == 2)
        if cpGen.params[1] > 0
          objPar1p[i] = 1/sqrt(cpGen.params[1]*fData.baseMVA);
          objPar2p[i] = sqrt(cpGen.params[1]*fData.baseMVA);
        else
          objPar1p[i] = 1;
          objPar2p[i] = 1;
        end
        # linear or constant
        @constraint(mp,tp[i] == objPar1p[i]*(cpGen.params[1]*fData.baseMVA*sps[i]) + cpGen.params[2]);
        @constraint(mp,sps[i] == sp[i]*objPar2p[i]);
      elseif (cpGen.n == 1)
        @constraint(mp,tp[i] == cpGen.params[1]);
      end
      objExpr += tp[i];
    end
    if !(fData.cq == Dict())
      cqGen = fData.cq[i];
      if cqGen.model == 1
        if cqGen.n > 4
          xq[i] = @variable(mp,[1:(div(cqGen.n,2)-1)],Bin);
          zq[i] = @variable(mp);
          @constraint(mp,sq[i] <= cqGen.params[cqGen.n-1]);
          @constraint(mp,sq[i] >= cqGen.params[1]);
          @constraint(mp,sum(xq[i][j] for j in 1:(div(cqGen.n,2)-1)) == 1);
          @constraint(mp,[j in 1:(div(cqGen.n,2)-1)],sq[i] <= (1-xq[i][j])*(cqGen.params[cqGen.n-1]-cqGen.params[1]+10)+cqGen.qarams[(j+1)*2]);
          @constraint(mp,[j in 1:(div(cqGen.n,2)-1)],sq[i] >= (xp[i][j]-1)*(cqGen.params[cqGen.n-1]-cqGen.params[1]+10)+cqGen.params[j*2]);
          @constraint(mp,[j in 1:(div(cqGen.n,2)-1)],zq[i] >= sq[i]*(cqGen.params[(j+1)*2] - cqGen.params[j*2])/
            (cqGen.params[(j+1)*2-1] - cqGen.params[j*2-1])+(cqGen.params[j*2]-(cqGen.params[(j+1)*2] - cqGen.params[j*2])/
            (cqGen.params[(j+1)*2-1] - cqGen.params[j*2-1])*cqGen.params[j*2-1]));
          objExpr += zq[j];
        elseif cqGen.n == 4
          @constraint(mp,sq[i] <= cqGen.params[cqGen.n-1]);
          @constraint(mp,sq[i] >= cqGen.params[1]);
          objExpr += (cqGen.params[4] - cqGen.params[2])/(cqGen.params[3] - cqGen.params[1])*sq[i] + cqGen.params[2];
        else
          error("Not enough parameters");
        end
      else
        if (cqGen.n == 3)&(cqGen.params[1] > 0)
          # quadratic
          αi = sqrt(cqGen.params[1])*fData.baseMVA;
          βi = cqGen.params[2]/(2*sqrt(cqGen.params[1]));
          f1i = -1/4 + cqGen.params[2]^2/(4*cqGen.params[1]) - cqGen.params[3];
          f2i = 1/4 + cqGen.params[2]^2/(4*cqGen.params[1]) - cqGen.params[3];
          @constraint(mp,norm([αi*sq[i] + βi, tq[i] + f1i]) <= tq[i] + f2i);
          @constraint(mp,tq[i] + f2i >= 0);
        elseif (cqGen.n == 3)&(cqGen.params[1] == 0)
          # quadratic but the quadratic coefficient is 0, so linear
          if cqGen.params[2] > 0
            objPar1q[i] = 1/sqrt(cqGen.params[2]*fData.baseMVA);
            objPar2q[i] = sqrt(cqGen.params[2]*fData.baseMVA);
          else
            objPar1q[i] = 1;
            objPar2q[i] = 1;
          end
          @constraint(mp,tq[i] == objPar1q[i]*(cqGen.params[2]*fData.baseMVA*sqs[i]) + cqGen.params[3]);
          @constraint(mp,sqs[i] == sq[i]*objPar2q[i]);
        elseif (cqGen.n == 2)
          # linear or constant
          if cqGen.params[1] > 0
            objPar1q[i] = 1/sqrt(cqGen.params[1]*fData.baseMVA);
            objPar2q[i] = sqrt(cqGen.params[1]*fData.baseMVA);
          else
            objPar1q[i] = 1;
            objPar2q[i] = 1;
          end
          @constraint(mp,tq[i] == objPar1q[i]*(cqGen.params[1]*sq[i]*fData.baseMVA) + cqGen.params[2]);
          @constraint(mp, sqs[i] == sq[i]*objPar2q[i]);
        elseif (cqGen.n == 1)
          @constraint(mp,tq[i] == cqGen.params[1]);
        end
        objExpr += tq[i];
      end
    end
  end
  @objective(mp,Min,objExpr);

  @variable(mp,p[k in fData.brList]);
  @variable(mp,q[k in fData.brList]);
  @variable(mp,ps[k in fData.brList]);
  @variable(mp,qs[k in fData.brList]);
  @variable(mp,vmin[i] <= v[i in fData.IDList] <= vmax[i]);
  @variable(mp,vhat[i in fData.IDList]);
  @variable(mp,θ[i in fData.IDList]);

  @variable(mp,opplus[i in fData.IDList] >= 0);
  @variable(mp,opminus[i in fData.IDList] >= 0);
  @variable(mp,oqplus[i in fData.IDList] >= 0);
  @variable(mp,oqminus[i in fData.IDList] >= 0);

  @variable(mp,tAux2[k in fData.brList]);
  @variable(mp,tAux3[k in fData.brList]);
  @variable(mp,tAux4[k in fData.brList] >= 0);
  @variable(mp,tAux5[i in fData.IDList]);
  @variable(mp,tAux6[i in fData.IDList] >= 0);

  @variable(mp,vv[k in fData.brList]);
  @variable(mp,cos(θu[k]) <= cs[k in fData.brList] <= 1);
  @variable(mp,-sin(θu[k]) <= ss[k in fData.brList] <= sin(θu[k]));
  # add the strengthened bounds on cs
  csmax = Dict();
  csmin = Dict();
  ssmax = Dict();
  ssmin = Dict();
  for k in fData.brList
    if (θDmax[k] >= 0)&(θDmin[k] >= 0)
      csmax[k] = cos(θDmin[k]);
      csmin[k] = cos(θDmax[k]);
      @constraint(mp, cs[k] <= csmax[k]);
      @constraint(mp, cs[k] >= csmin[k]);
    elseif (θDmax[k] < 0)&(θDmin[k] < 0)
      csmax[k] = cos(θDmax[k]);
      csmin[k] = cos(θDmin[k]);
      @constraint(mp, cs[k] <= csmax[k]);
      @constraint(mp, cs[k] >= csmin[k]);
    else
      csmax[k] = 1;
      csmin[k] = min(cos(θDmax[k]),cos(θDmin[k]));
      @constraint(mp, cs[k] <= csmax[k]);
      @constraint(mp, cs[k] >= csmin[k]);
    end
    @constraint(mp, cs[k] >= (cos(θDmax[k]) - cos(θDmin[k]))/(θDmax[k] - θDmin[k])*((θ[k[1]] - fData.σ[k]) - θ[k[2]] - θDmin[k]) + cos(θDmin[k]));
  end
  # add the strengthened bounds on ss
  for k in fData.brList
    ssmax[k] = sin(θDmax[k]);
    ssmin[k] = sin(θDmin[k]);
    @constraint(mp, ss[k] <= ssmax[k]);
    @constraint(mp, ss[k] >= ssmin[k]);
  end

  @variable(mp,wc[k in fData.brList]);
  @variable(mp,ws[k in fData.brList]);
  @constraint(mp,wcEquality[k in fData.brList;k[1] < k[2]], wc[k] == wc[(k[2],k[1],k[3])]);
  @constraint(mp,wsEquality[k in fData.brList;k[1] < k[2]], ws[k] == -ws[(k[2],k[1],k[3])]);
  @constraint(mp,vvEquality[k in fData.brList;k[1] < k[2]], vv[k] == vv[(k[2],k[1],k[3])]);
  @constraint(mp,csEquality[k in fData.brList;k[1] < k[2]], cs[k] == cs[(k[2],k[1],k[3])]);
  @constraint(mp,ssEquality[k in fData.brList;k[1] < k[2]], ss[k] == -ss[(k[2],k[1],k[3])]);

  con1Par = Dict();
  con2Par = Dict();
  for k in fData.brList
    con1Par[k] = 1/sqrt(abs(fData.b[k]));
    con2Par[k] = sqrt(abs(fData.b[k]));
  end

  @constraint(mp,lineConstrP[k in fData.brList],ps[k] == con1Par[k]*(fData.g[k]*vhat[k[1]]/(fData.τ1[k]^2) - fData.g[k]*wc[k] - fData.b[k]*ws[k]));
  @constraint(mp,pScale[k in fData.brList],p[k] == con2Par[k]*ps[k]);
  @constraint(mp,lineConstrQ[k in fData.brList],qs[k] == con1Par[k]*((-fData.b[k] - fData.bc[k]/2)*vhat[k[1]]/(fData.τ1[k]^2) + fData.b[k]*wc[k] - fData.g[k]*ws[k]));
  @constraint(mp,qScale[k in fData.brList],q[k] == con2Par[k]*qs[k]);
  socList1 = Dict();
  socList2 = Dict();
  socList3 = Dict();
  socList4 = Dict();
  socList5 = Dict();
  for k in fData.brList
    if fData.rateA[k] < Inf
      socList1[k] = [p[k],q[k]];
    end
    socList2[k] = [vv[k],vhat[k[1]]/((fData.τ1[k]^2)*sqrt(2)),vhat[k[2]]/((fData.τ2[k]^2)*sqrt(2))];
    socList3[k] = [tAux2[k],tAux3[k]];
    socList5[k] = [wc[k],ws[k],vhat[k[1]]/((fData.τ1[k]^2)*sqrt(2)),vhat[k[2]]/((fData.τ2[k]^2)*sqrt(2))];
  end
  for i in fData.IDList
    socList4[i] = [v[i],tAux5[i]];
  end

  @constraint(mp,socConstraint1[k in fData.brList; fData.rateA[k] < Inf],norm(socList1[k]) <= fData.rateA[k]);
  @constraint(mp,socConstraint2[k in fData.brList], norm(socList2[k]) <= (vhat[k[1]]/(fData.τ1[k]^2) + vhat[k[2]]/(fData.τ2[k]^2))/sqrt(2));
  @constraint(mp,socConstraint3[k in fData.brList],norm(socList3[k]) <= tAux4[k]);
  @constraint(mp,socConstraint4[i in fData.IDList],norm(socList4[i]) <= tAux6[i]);
  @constraint(mp,socConstraint5[k in fData.brList],norm(socList5[k]) <= (vhat[k[1]]/(fData.τ1[k]^2) + vhat[k[2]]/(fData.τ2[k]^2))/sqrt(2));

  @constraint(mp,auxConstr2[k in fData.brList],sqrt((1-cos(θu[k]))/(θu[k])^2)*((θ[k[1]] - fData.σ[k]) - θ[k[2]]) == tAux2[k]);
  @constraint(mp,auxConstr3[k in fData.brList],cs[k] - 3/4 == tAux3[k]);
  @constraint(mp,auxConstr4[k in fData.brList],5/4 - cs[k] == tAux4[k]);
  @constraint(mp,sinMock1[k in fData.brList],
              ss[k] <= cos(θu[k]/2)*((θ[k[1]] - fData.σ[k]) - θ[k[2]] - θu[k]/2) + sin(θu[k]/2));
  @constraint(mp,sinMock2[k in fData.brList],
              ss[k] >= cos(θu[k]/2)*((θ[k[1]] - fData.σ[k]) - θ[k[2]] + θu[k]/2) - sin(θu[k]/2));
  @constraint(mp,angleDiff1[k in fData.brList],(θ[k[1]] - fData.σ[k]) - θ[k[2]] <= θDmax[k]);
  @constraint(mp,angleDiff2[k in fData.brList],(θ[k[1]] - fData.σ[k]) - θ[k[2]] >= θDmin[k]);
  @constraint(mp,auxConstr5[i in fData.IDList],tAux5[i] == vhat[i] - 1/4);
  @constraint(mp,auxConstr6[i in fData.IDList],tAux6[i] == vhat[i] + 1/4);
  @constraint(mp,v2Mock2[i in fData.IDList], vhat[i] - (vmax[i] + vmin[i])*v[i] <= -vmax[i]*vmin[i]);
  @constraint(mp,vvMock1[k in fData.brList],
              vv[k] >= vmin[k[1]]*v[k[2]]/(fData.τ1[k]*fData.τ2[k]) + vmin[k[2]]*v[k[1]]/(fData.τ1[k]*fData.τ2[k]) - vmin[k[1]]*vmin[k[2]]/(fData.τ1[k]*fData.τ2[k]));
  @constraint(mp,vvMock2[k in fData.brList],
              vv[k] >= vmax[k[1]]*v[k[2]]/(fData.τ1[k]*fData.τ2[k]) + vmax[k[2]]*v[k[1]]/(fData.τ1[k]*fData.τ2[k]) - vmax[k[1]]*vmax[k[2]]/(fData.τ1[k]*fData.τ2[k]));
  @constraint(mp,vvMock3[k in fData.brList],
              vv[k] <= vmin[k[1]]*v[k[2]]/(fData.τ1[k]*fData.τ2[k]) + vmax[k[2]]*v[k[1]]/(fData.τ1[k]*fData.τ2[k]) - vmin[k[1]]*vmax[k[2]]/(fData.τ1[k]*fData.τ2[k]));
  @constraint(mp,vvMock4[k in fData.brList],
              vv[k] <= vmax[k[1]]*v[k[2]]/(fData.τ1[k]*fData.τ2[k]) + vmin[k[2]]*v[k[1]]/(fData.τ1[k]*fData.τ2[k]) - vmax[k[1]]*vmin[k[2]]/(fData.τ1[k]*fData.τ2[k]));

  @constraint(mp,wcMock1[k in fData.brList],
              wc[k] >= vmin[k[1]]*vmin[k[2]]/(fData.τ1[k]*fData.τ2[k])*cs[k] + csmin[k]*vv[k] - vmin[k[1]]*vmin[k[2]]/(fData.τ1[k]*fData.τ2[k])*csmin[k]);
  @constraint(mp,wcMock2[k in fData.brList],
              wc[k] >= vmax[k[1]]*vmax[k[2]]/(fData.τ1[k]*fData.τ2[k])*cs[k] + csmax[k]*vv[k] - vmax[k[1]]*vmax[k[2]]/(fData.τ1[k]*fData.τ2[k])*csmax[k]);
  @constraint(mp,wcMock3[k in fData.brList],
              wc[k] <= vmin[k[1]]*vmin[k[2]]/(fData.τ1[k]*fData.τ2[k])*cs[k] + vv[k]*csmax[k] - vmin[k[1]]*vmin[k[2]]/(fData.τ1[k]*fData.τ2[k])*csmax[k]);
  @constraint(mp,wcMock4[k in fData.brList],
              wc[k] <= vmax[k[1]]*vmax[k[2]]/(fData.τ1[k]*fData.τ2[k])*cs[k] + vv[k]*csmin[k] - vmax[k[1]]*vmax[k[2]]/(fData.τ1[k]*fData.τ2[k])*csmin[k]);

  @constraint(mp,wsMock1[k in fData.brList],
              ws[k] >= vmin[k[1]]*vmin[k[2]]/(fData.τ1[k]*fData.τ2[k])*ss[k] + ssmin[k]*vv[k] - vmin[k[1]]*vmin[k[2]]/(fData.τ1[k]*fData.τ2[k])*ssmin[k]);
  @constraint(mp,wsMock2[k in fData.brList],
              ws[k] >= vmax[k[1]]*vmax[k[2]]/(fData.τ1[k]*fData.τ2[k])*ss[k] + ssmax[k]*vv[k] - vmax[k[1]]*vmax[k[2]]/(fData.τ1[k]*fData.τ2[k])*ssmax[k]);
  @constraint(mp,wsMock3[k in fData.brList],
              ws[k] <= vmin[k[1]]*vmin[k[2]]/(fData.τ1[k]*fData.τ2[k])*ss[k] + vv[k]*ssmax[k] - vmin[k[1]]*vmin[k[2]]/(fData.τ1[k]*fData.τ2[k])*ssmax[k]);
  @constraint(mp,wsMock4[k in fData.brList],
              ws[k] <= vmax[k[1]]*vmax[k[2]]/(fData.τ1[k]*fData.τ2[k])*ss[k] + vv[k]*ssmin[k] - vmax[k[1]]*vmax[k[2]]/(fData.τ1[k]*fData.τ2[k])*ssmin[k]);

  @constraint(mp,tanConstr1[k in fData.brList],
                ws[k] - tan(θDmax[k])*wc[k] <= 0);
  @constraint(mp,tanConstr2[k in fData.brList],
                ws[k] - tan(θDmin[k])*wc[k] >= 0);

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

  @constraint(mp,lncConstr1[k in fData.brList],
              vδ[k[1]]*vδ[k[2]]/(fData.τ1[k]*fData.τ2[k])*(wc[k]*cos(θϕ[k]) + ws[k]*sin(θϕ[k])) -
              vmax[k[2]]/fData.τ2[k]*cos(θδ[k])*vδ[k[2]]/fData.τ2[k]*vhat[k[1]]/(fData.τ1[k]^2) -
              vmax[k[1]]/fData.τ1[k]*cos(θδ[k])*vδ[k[1]]/fData.τ1[k]*vhat[k[2]]/(fData.τ2[k]^2) >=
              vmax[k[1]]*vmax[k[2]]/(fData.τ1[k]*fData.τ2[k])*cos(θδ[k])*(vmin[k[1]]*vmin[k[2]]/(fData.τ1[k]*fData.τ2[k]) -
              vmax[k[1]]*vmax[k[2]]/(fData.τ1[k]*fData.τ2[k])));
  @constraint(mp,lncConstr2[k in fData.brList],
              vδ[k[1]]*vδ[k[2]]/(fData.τ1[k]*fData.τ2[k])*(wc[k]*cos(θϕ[k]) + ws[k]*sin(θϕ[k])) -
              vmin[k[2]]/fData.τ2[k]*cos(θδ[k])*vδ[k[2]]/fData.τ2[k]*vhat[k[1]]/(fData.τ1[k]^2) -
              vmin[k[1]]/fData.τ1[k]*cos(θδ[k])*vδ[k[1]]/fData.τ1[k]*vhat[k[2]]/(fData.τ2[k]^2) >=
              -vmin[k[1]]*vmin[k[2]]/(fData.τ1[k]*fData.τ2[k])*cos(θδ[k])*(vmin[k[1]]*vmin[k[2]]/(fData.τ1[k]*fData.τ2[k]) -
              vmax[k[1]]*vmax[k[2]]/(fData.τ1[k]*fData.τ2[k])));

  conWPpar1 = Dict();
  conWQpar1 = Dict();
  conWPpar2 = Dict();
  conWQpar2 = Dict();
  for i in fData.IDList
    if (abs(uData[i].WPmax) <= 1e-3)&(abs(uData[i].WPmax) < 0)
      conWPpar1[i] = 1/sqrt(abs(uData[i].WPmax));
      conWPpar2[i] = sqrt(abs(uData[i].WPmax));
    else
      conWPpar1[i] = 1;
      conWPpar2[i] = 1;
    end
    if (abs(uData[i].WQmax) <= 1e-3)&(abs(uData[i].WQmax) < 0)
      conWQpar1[i] = 1/sqrt(abs(uData[i].WQmax));
      conWQpar2[i] = sqrt(abs(uData[i].WQmax));
    else
      conWQpar1[i] = 1;
      conWQpar2[i] = 1;
    end
  end
  if βSwitch == 1
    @constraint(mp,opConstr1[i in fData.IDList],opminus[i] <= conWQpar1[i]*uData[i].WPmax);
    @constraint(mp,oqConstr1[i in fData.IDList],oqminus[i] <= conWQpar1[i]*uData[i].WQmax);
    @constraint(mp,opConstr2[i in fData.IDList],opplus[i] <= conWQpar1[i]*(uData[i].WPmax + uData[i].RESP0));
    @constraint(mp,oqConstr2[i in fData.IDList],oqplus[i] <= conWQpar1[i]*(uData[i].WQmax + uData[i].RESQ0));
  else
    @constraint(mp,opConstr[i in fData.IDList],opplus[i] + opminus[i] <= 0);
    @constraint(mp,oqConstr[i in fData.IDList],oqplus[i] + oqminus[i] <= 0);
  end

  @constraint(mp,totalP[i in fData.IDList], sum(p[k] for k in branchDict1[i]) + vhat[i]*fData.gs[i]
      + opplus[i]*conWPpar2[i] - opminus[i]*conWPpar2[i] == sphatsum[i] + uData[i].GP0);
  @constraint(mp,totalQ[i in fData.IDList], sum(q[k] for k in branchDict1[i]) - vhat[i]*fData.bs[i]
      + oqplus[i]*conWQpar2[i] - oqminus[i]*conWQpar2[i] == sqhatsum[i] + uData[i].GQ0);

  return mp;
end

function ACOPF_C_QP(fData,uData,vmax,vmin,θDmax,θDmin,βSwitch = 1)
  θu = Dict();
  for k in fData.brList
    θu[k] = max(abs(θDmax[k]),abs(θDmin[k]));
  end
  # scaling issue exists for this primal problem.
  mp = Model(solver = IpoptSolver(print_level = 0,max_iter = 10000,linear_solver = "ma27"));
  # mp = Model(solver = MosekSolver(MSK_DPAR_MIO_TOL_REL_GAP = 1e-7,MSK_DPAR_INTPNT_TOL_PFEAS = 1e-9,MSK_DPAR_MIO_TOL_FEAS = 1e-9,MSK_IPAR_LOG = 0));

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

  @variable(mp,p[k in fData.brList]);
  @variable(mp,q[k in fData.brList]);
  @variable(mp,ps[k in fData.brList]);
  @variable(mp,qs[k in fData.brList]);
  @variable(mp,vmin[i] <= v[i in fData.IDList] <= vmax[i]);
  @variable(mp,vhat[i in fData.IDList]);
  @variable(mp,θ[i in fData.IDList]);

  @variable(mp,opplus[i in fData.IDList] >= 0);
  @variable(mp,opminus[i in fData.IDList] >= 0);
  @variable(mp,oqplus[i in fData.IDList] >= 0);
  @variable(mp,oqminus[i in fData.IDList] >= 0);

  @variable(mp,vv[k in fData.brList]);
  @variable(mp,cos(θu[k]) <= cs[k in fData.brList] <= 1);
  @variable(mp,-sin(θu[k]) <= ss[k in fData.brList] <= sin(θu[k]));
  # add the strengthened bounds on cs
  csmax = Dict();
  csmin = Dict();
  ssmax = Dict();
  ssmin = Dict();
  for k in fData.brList
    if (θDmax[k] >= 0)&(θDmin[k] >= 0)
      csmax[k] = cos(θDmin[k]);
      csmin[k] = cos(θDmax[k]);
      @constraint(mp, cs[k] <= csmax[k]);
      @constraint(mp, cs[k] >= csmin[k]);
    elseif (θDmax[k] < 0)&(θDmin[k] < 0)
      csmax[k] = cos(θDmax[k]);
      csmin[k] = cos(θDmin[k]);
      @constraint(mp, cs[k] <= csmax[k]);
      @constraint(mp, cs[k] >= csmin[k]);
    else
      csmax[k] = 1;
      csmin[k] = min(cos(θDmax[k]),cos(θDmin[k]));
      @constraint(mp, cs[k] <= csmax[k]);
      @constraint(mp, cs[k] >= csmin[k]);
    end
    @constraint(mp, cs[k] >= (cos(θDmax[k]) - cos(θDmin[k]))/(θDmax[k] - θDmin[k])*((θ[k[1]] - fData.σ[k]) - θ[k[2]] - θDmin[k]) + cos(θDmin[k]));
  end
  # add the strengthened bounds on ss
  for k in fData.brList
    ssmax[k] = sin(θDmax[k]);
    ssmin[k] = sin(θDmin[k]);
    @constraint(mp, ss[k] <= ssmax[k]);
    @constraint(mp, ss[k] >= ssmin[k]);
  end

  @variable(mp,wc[k in fData.brList]);
  @variable(mp,ws[k in fData.brList]);
  @constraint(mp,wcEquality[k in fData.brList;k[1] < k[2]], wc[k] == wc[(k[2],k[1],k[3])]);
  @constraint(mp,wsEquality[k in fData.brList;k[1] < k[2]], ws[k] == -ws[(k[2],k[1],k[3])]);
  @constraint(mp,vvEquality[k in fData.brList;k[1] < k[2]], vv[k] == vv[(k[2],k[1],k[3])]);
  @constraint(mp,csEquality[k in fData.brList;k[1] < k[2]], cs[k] == cs[(k[2],k[1],k[3])]);
  @constraint(mp,ssEquality[k in fData.brList;k[1] < k[2]], ss[k] == -ss[(k[2],k[1],k[3])]);

  con1Par = Dict();
  con2Par = Dict();
  for k in fData.brList
    con1Par[k] = 1/sqrt(abs(fData.b[k]));
    con2Par[k] = sqrt(abs(fData.b[k]));
  end

  @constraint(mp,lineConstrP[k in fData.brList],ps[k] == con1Par[k]*(fData.g[k]*vhat[k[1]]/(fData.τ1[k]^2) - fData.g[k]*wc[k] - fData.b[k]*ws[k]));
  @constraint(mp,pScale[k in fData.brList],p[k] == con2Par[k]*ps[k]);
  @constraint(mp,lineConstrQ[k in fData.brList],qs[k] == con1Par[k]*((-fData.b[k] - fData.bc[k]/2)*vhat[k[1]]/(fData.τ1[k]^2) + fData.b[k]*wc[k] - fData.g[k]*ws[k]));
  @constraint(mp,qScale[k in fData.brList],q[k] == con2Par[k]*qs[k]);

  @constraint(mp,socConstraint1[k in fData.brList; fData.rateA[k] < Inf],p[k]^2 + q[k]^2 <= fData.rateA[k]^2);
  @NLconstraint(mp,socConstraint2[k in fData.brList],vv[k]^2 <= vhat[k[1]]/fData.τ1[k]^2 * vhat[k[2]]/fData.τ2[k]^2);
  @NLconstraint(mp,socConstraint3[k in fData.brList],cs[k] + (1-cos(θu[k]))/(θu[k])^2*((θ[k[1]] - fData.σ[k]) - θ[k[2]])^2 <= 1);
  @NLconstraint(mp,socConstraint4[i in fData.IDList],v[i]^2 <= vhat[i]);
  @NLconstraint(mp,socConstraint5[k in fData.brList],wc[k]^2 + ws[k]^2 <= vhat[k[1]]/(fData.τ1[k]^2)*vhat[k[2]]/(fData.τ2[k]^2));

  @constraint(mp,sinMock1[k in fData.brList],
              ss[k] <= cos(θu[k]/2)*((θ[k[1]] - fData.σ[k]) - θ[k[2]] - θu[k]/2) + sin(θu[k]/2));
  @constraint(mp,sinMock2[k in fData.brList],
              ss[k] >= cos(θu[k]/2)*((θ[k[1]] - fData.σ[k]) - θ[k[2]] + θu[k]/2) - sin(θu[k]/2));
  @constraint(mp,angleDiff1[k in fData.brList],(θ[k[1]] - fData.σ[k]) - θ[k[2]] <= θDmax[k]);
  @constraint(mp,angleDiff2[k in fData.brList],(θ[k[1]] - fData.σ[k]) - θ[k[2]] >= θDmin[k]);
  @constraint(mp,v2Mock2[i in fData.IDList], vhat[i] - (vmax[i] + vmin[i])*v[i] <= -vmax[i]*vmin[i]);
  @constraint(mp,vvMock1[k in fData.brList],
              vv[k] >= vmin[k[1]]*v[k[2]]/(fData.τ1[k]*fData.τ2[k]) + vmin[k[2]]*v[k[1]]/(fData.τ1[k]*fData.τ2[k]) - vmin[k[1]]*vmin[k[2]]/(fData.τ1[k]*fData.τ2[k]));
  @constraint(mp,vvMock2[k in fData.brList],
              vv[k] >= vmax[k[1]]*v[k[2]]/(fData.τ1[k]*fData.τ2[k]) + vmax[k[2]]*v[k[1]]/(fData.τ1[k]*fData.τ2[k]) - vmax[k[1]]*vmax[k[2]]/(fData.τ1[k]*fData.τ2[k]));
  @constraint(mp,vvMock3[k in fData.brList],
              vv[k] <= vmin[k[1]]*v[k[2]]/(fData.τ1[k]*fData.τ2[k]) + vmax[k[2]]*v[k[1]]/(fData.τ1[k]*fData.τ2[k]) - vmin[k[1]]*vmax[k[2]]/(fData.τ1[k]*fData.τ2[k]));
  @constraint(mp,vvMock4[k in fData.brList],
              vv[k] <= vmax[k[1]]*v[k[2]]/(fData.τ1[k]*fData.τ2[k]) + vmin[k[2]]*v[k[1]]/(fData.τ1[k]*fData.τ2[k]) - vmax[k[1]]*vmin[k[2]]/(fData.τ1[k]*fData.τ2[k]));

  @constraint(mp,wcMock1[k in fData.brList],
              wc[k] >= vmin[k[1]]*vmin[k[2]]/(fData.τ1[k]*fData.τ2[k])*cs[k] + csmin[k]*vv[k] - vmin[k[1]]*vmin[k[2]]/(fData.τ1[k]*fData.τ2[k])*csmin[k]);
  @constraint(mp,wcMock2[k in fData.brList],
              wc[k] >= vmax[k[1]]*vmax[k[2]]/(fData.τ1[k]*fData.τ2[k])*cs[k] + csmax[k]*vv[k] - vmax[k[1]]*vmax[k[2]]/(fData.τ1[k]*fData.τ2[k])*csmax[k]);
  @constraint(mp,wcMock3[k in fData.brList],
              wc[k] <= vmin[k[1]]*vmin[k[2]]/(fData.τ1[k]*fData.τ2[k])*cs[k] + vv[k]*csmax[k] - vmin[k[1]]*vmin[k[2]]/(fData.τ1[k]*fData.τ2[k])*csmax[k]);
  @constraint(mp,wcMock4[k in fData.brList],
              wc[k] <= vmax[k[1]]*vmax[k[2]]/(fData.τ1[k]*fData.τ2[k])*cs[k] + vv[k]*csmin[k] - vmax[k[1]]*vmax[k[2]]/(fData.τ1[k]*fData.τ2[k])*csmin[k]);

  @constraint(mp,wsMock1[k in fData.brList],
              ws[k] >= vmin[k[1]]*vmin[k[2]]/(fData.τ1[k]*fData.τ2[k])*ss[k] + ssmin[k]*vv[k] - vmin[k[1]]*vmin[k[2]]/(fData.τ1[k]*fData.τ2[k])*ssmin[k]);
  @constraint(mp,wsMock2[k in fData.brList],
              ws[k] >= vmax[k[1]]*vmax[k[2]]/(fData.τ1[k]*fData.τ2[k])*ss[k] + ssmax[k]*vv[k] - vmax[k[1]]*vmax[k[2]]/(fData.τ1[k]*fData.τ2[k])*ssmax[k]);
  @constraint(mp,wsMock3[k in fData.brList],
              ws[k] <= vmin[k[1]]*vmin[k[2]]/(fData.τ1[k]*fData.τ2[k])*ss[k] + vv[k]*ssmax[k] - vmin[k[1]]*vmin[k[2]]/(fData.τ1[k]*fData.τ2[k])*ssmax[k]);
  @constraint(mp,wsMock4[k in fData.brList],
              ws[k] <= vmax[k[1]]*vmax[k[2]]/(fData.τ1[k]*fData.τ2[k])*ss[k] + vv[k]*ssmin[k] - vmax[k[1]]*vmax[k[2]]/(fData.τ1[k]*fData.τ2[k])*ssmin[k]);

  @constraint(mp,tanConstr1[k in fData.brList],
                ws[k] - tan(θDmax[k])*wc[k] <= 0);
  @constraint(mp,tanConstr2[k in fData.brList],
                ws[k] - tan(θDmin[k])*wc[k] >= 0);

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

  @constraint(mp,lncConstr1[k in fData.brList],
              vδ[k[1]]*vδ[k[2]]/(fData.τ1[k]*fData.τ2[k])*(wc[k]*cos(θϕ[k]) + ws[k]*sin(θϕ[k])) -
              vmax[k[2]]/fData.τ2[k]*cos(θδ[k])*vδ[k[2]]/fData.τ2[k]*vhat[k[1]]/(fData.τ1[k]^2) -
              vmax[k[1]]/fData.τ1[k]*cos(θδ[k])*vδ[k[1]]/fData.τ1[k]*vhat[k[2]]/(fData.τ2[k]^2) >=
              vmax[k[1]]*vmax[k[2]]/(fData.τ1[k]*fData.τ2[k])*cos(θδ[k])*(vmin[k[1]]*vmin[k[2]]/(fData.τ1[k]*fData.τ2[k]) -
              vmax[k[1]]*vmax[k[2]]/(fData.τ1[k]*fData.τ2[k])));
  @constraint(mp,lncConstr2[k in fData.brList],
              vδ[k[1]]*vδ[k[2]]/(fData.τ1[k]*fData.τ2[k])*(wc[k]*cos(θϕ[k]) + ws[k]*sin(θϕ[k])) -
              vmin[k[2]]/fData.τ2[k]*cos(θδ[k])*vδ[k[2]]/fData.τ2[k]*vhat[k[1]]/(fData.τ1[k]^2) -
              vmin[k[1]]/fData.τ1[k]*cos(θδ[k])*vδ[k[1]]/fData.τ1[k]*vhat[k[2]]/(fData.τ2[k]^2) >=
              -vmin[k[1]]*vmin[k[2]]/(fData.τ1[k]*fData.τ2[k])*cos(θδ[k])*(vmin[k[1]]*vmin[k[2]]/(fData.τ1[k]*fData.τ2[k]) -
              vmax[k[1]]*vmax[k[2]]/(fData.τ1[k]*fData.τ2[k])));

  conWPpar1 = Dict();
  conWQpar1 = Dict();
  conWPpar2 = Dict();
  conWQpar2 = Dict();
  for i in fData.IDList
    if (abs(uData[i].WPmax) <= 1e-3)&(abs(uData[i].WPmax) < 0)
      conWPpar1[i] = 1/sqrt(abs(uData[i].WPmax));
      conWPpar2[i] = sqrt(abs(uData[i].WPmax));
    else
      conWPpar1[i] = 1;
      conWPpar2[i] = 1;
    end
    if (abs(uData[i].WQmax) <= 1e-3)&(abs(uData[i].WQmax) < 0)
      conWQpar1[i] = 1/sqrt(abs(uData[i].WQmax));
      conWQpar2[i] = sqrt(abs(uData[i].WQmax));
    else
      conWQpar1[i] = 1;
      conWQpar2[i] = 1;
    end
  end
  if βSwitch == 1
    @constraint(mp,opConstr1[i in fData.IDList],opminus[i] <= conWQpar1[i]*uData[i].WPmax);
    @constraint(mp,oqConstr1[i in fData.IDList],oqminus[i] <= conWQpar1[i]*uData[i].WQmax);
    @constraint(mp,opConstr2[i in fData.IDList],opplus[i] <= conWQpar1[i]*(uData[i].WPmax + uData[i].RESP0));
    @constraint(mp,oqConstr2[i in fData.IDList],oqplus[i] <= conWQpar1[i]*(uData[i].WQmax + uData[i].RESQ0));
  else
    @constraint(mp,opConstr[i in fData.IDList],opplus[i] + opminus[i] <= 0);
    @constraint(mp,oqConstr[i in fData.IDList],oqplus[i] + oqminus[i] <= 0);
  end

  @constraint(mp,totalP[i in fData.IDList], sum(p[k] for k in branchDict1[i]) + vhat[i]*fData.gs[i]
      + opplus[i]*conWPpar2[i] - opminus[i]*conWPpar2[i] == sphatsum[i] + uData[i].GP0);
  @constraint(mp,totalQ[i in fData.IDList], sum(q[k] for k in branchDict1[i]) - vhat[i]*fData.bs[i]
      + oqplus[i]*conWQpar2[i] - oqminus[i]*conWQpar2[i] == sqhatsum[i] + uData[i].GQ0);

  return mp;
end
