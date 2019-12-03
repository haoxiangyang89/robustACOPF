# append the entire scenario to the first stage
function appendScen(fData,uData,Ω,yp,yq,vmax,vmin,θDmax,θDmin,cutList,numericOpt = 0,QP_opt = 0)
  θu = Dict();
  for k in fData.brList
    θu[k] = max(abs(θDmax[k]),abs(θDmin[k]));
  end

  # scaling issue exists for this primal problem.
  if QP_opt == 0
    mp = Model(solver = GurobiSolver(NumericFocus = numericOpt,OutputFlag = 0, Threads = 20));
  else
    mp = Model(solver = IpoptSolver(print_level = 0,linear_solver = "ma27"));
  end

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

  # construct the objective function
  @expression(mp,objExpr,0);
  xp = Dict();
  xq = Dict();
  zp = Dict();
  zq = Dict();
  @variable(mp,tp[i in fData.genIDList]);
  @variable(mp,tq[i in fData.genIDList]);
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
      if QP_opt == 0
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
          @constraint(mp,tp[i] == cpGen.params[2]*sp[i]*fData.baseMVA + cpGen.params[3]);
        elseif (cpGen.n == 2)
          # linear or constant
          @constraint(mp,tp[i] == cpGen.params[1]*sp[i]*fData.baseMVA + cpGen.params[2]);
        elseif (cpGen.n == 1)
          @constraint(mp,tp[i] == cpGen.params[1]);
        end
        objExpr += tp[i];
      else
        objExpr += cpGen.params[1]*(sp[i]*fData.baseMVA)^2 + cpGen.params[2]*sp[i]*fData.baseMVA + cpGen.params[3];
      end
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
        if QP_opt == 0
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
            @constraint(mp,tq[i] == cqGen.params[2]*sq[i]*fData.baseMVA + cqGen.params[3]);
          elseif (cqGen.n == 2)
            # linear or constant
            @constraint(mp,tq[i] == cqGen.params[1]*sq[i]*fData.baseMVA + cqGen.params[2]);
          elseif (cqGen.n == 1)
            @constraint(mp,tq[i] == cqGen.params[1]);
          end
          objExpr += tq[i];
        else
          objExpr += cqGen.params[1]*(sq[i]*fData.baseMVA)^2 + cqGen.params[2]*sq[i]*fData.baseMVA + cqGen.params[3];
        end
      end
    end
  end
  @objective(mp,Min,objExpr);

  @variable(mp,opplus[i in fData.IDList,ω in Ω] >= 0);
  @variable(mp,opminus[i in fData.IDList,ω in Ω] >= 0);
  @variable(mp,oqplus[i in fData.IDList,ω in Ω] >= 0);
  @variable(mp,oqminus[i in fData.IDList,ω in Ω] >= 0);

  @variable(mp,p[k in fData.brList,ω in Ω]);
  @variable(mp,q[k in fData.brList,ω in Ω]);
  @variable(mp,ps[k in fData.brList,ω in Ω]);
  @variable(mp,qs[k in fData.brList,ω in Ω]);
  @variable(mp,vmin[i] <= v[i in fData.IDList,ω in Ω] <= vmax[i]);
  @variable(mp,vhat[i in fData.IDList,ω in Ω]);
  @variable(mp,θ[i in fData.IDList,ω in Ω]);
  for i in fData.IDList
    if fData.bType[i] == 3
      for ω in Ω
        @constraint(mp,θ[i,ω] == 0);
      end
    end
  end

  @variable(mp,tAux2[k in fData.brList,ω in Ω]);
  @variable(mp,tAux3[k in fData.brList,ω in Ω]);
  @variable(mp,tAux4[k in fData.brList,ω in Ω] >= 0);
  @variable(mp,tAux5[i in fData.IDList,ω in Ω]);
  @variable(mp,tAux6[i in fData.IDList,ω in Ω] >= 0);

  @variable(mp,vv[k in fData.brList,ω in Ω]);
  @variable(mp,cos(θu[k]) <= cs[k in fData.brList,ω in Ω] <= 1);
  @variable(mp,-sin(θu[k]) <= ss[k in fData.brList,ω in Ω] <= sin(θu[k]));
  # add the strengthened bounds on cs
  csmax = Dict();
  csmin = Dict();
  ssmax = Dict();
  ssmin = Dict();
  for k in fData.brList
    if (θDmax[k] >= 0)&(θDmin[k] >= 0)
      csmax[k] = cos(θDmin[k]);
      csmin[k] = cos(θDmax[k]);
    elseif (θDmax[k] < 0)&(θDmin[k] < 0)
      csmax[k] = cos(θDmax[k]);
      csmin[k] = cos(θDmin[k]);
    else
      csmax[k] = 1;
      csmin[k] = min(cos(θDmax[k]),cos(θDmin[k]));
    end
  end
  @constraint(mp, csCon1[k in fData.brList,ω in Ω], cs[k,ω] <= csmax[k]);
  @constraint(mp, csCon2[k in fData.brList,ω in Ω], cs[k,ω] >= csmin[k]);
  @constraint(mp, csCon3[k in fData.brList,ω in Ω], cs[k,ω] >= (cos(θDmax[k]) - cos(θDmin[k]))/(θDmax[k] - θDmin[k])*((θ[k[1],ω] - fData.σ[k]) - θ[k[2],ω] - θDmin[k]) + cos(θDmin[k]));
  # add the strengthened bounds on ss
  for k in fData.brList
    ssmax[k] = sin(θDmax[k]);
    ssmin[k] = sin(θDmin[k]);
  end
  @constraint(mp, ssCon1[k in fData.brList,ω in Ω], ss[k,ω] <= ssmax[k]);
  @constraint(mp, ssCon2[k in fData.brList,ω in Ω], ss[k,ω] >= ssmin[k]);

  @variable(mp,wc[k in fData.brList,ω in Ω]);
  @variable(mp,ws[k in fData.brList,ω in Ω]);
  @constraint(mp,wcEquality[k in fData.brList,ω in Ω;k[1] < k[2]], wc[k,ω] == wc[(k[2],k[1],k[3]),ω]);
  @constraint(mp,wsEquality[k in fData.brList,ω in Ω;k[1] < k[2]], ws[k,ω] == -ws[(k[2],k[1],k[3]),ω]);
  @constraint(mp,vvEquality[k in fData.brList,ω in Ω;k[1] < k[2]], vv[k,ω] == vv[(k[2],k[1],k[3]),ω]);
  @constraint(mp,csEquality[k in fData.brList,ω in Ω;k[1] < k[2]], cs[k,ω] == cs[(k[2],k[1],k[3]),ω]);
  @constraint(mp,ssEquality[k in fData.brList,ω in Ω;k[1] < k[2]], ss[k,ω] == -ss[(k[2],k[1],k[3]),ω]);

  con1Par = Dict();
  con2Par = Dict();
  for k in fData.brList
    if abs(fData.b[k]) >= 1000
      con1Par[k] = 1/sqrt(abs(fData.b[k]));
      con2Par[k] = sqrt(abs(fData.b[k]));
    else
      con1Par[k] = 1;
      con2Par[k] = 1;
    end
  end

  @constraint(mp,lineConstrP[k in fData.brList,ω in Ω],ps[k,ω] == con1Par[k]*(fData.g[k]*vhat[k[1],ω]/(fData.τ1[k]^2) - fData.g[k]*wc[k,ω] - fData.b[k]*ws[k,ω]));
  @constraint(mp,pScale[k in fData.brList,ω in Ω],p[k,ω] == con2Par[k]*ps[k,ω]);
  @constraint(mp,lineConstrQ[k in fData.brList,ω in Ω],qs[k,ω] == con1Par[k]*((-fData.b[k] - fData.bc[k]/2)*vhat[k[1],ω]/(fData.τ1[k]^2) + fData.b[k]*wc[k,ω] - fData.g[k]*ws[k,ω]));
  @constraint(mp,qScale[k in fData.brList,ω in Ω],q[k,ω] == con2Par[k]*qs[k,ω]);
  if QP_opt == 0
    socList1 = Dict();
    socList2 = Dict();
    socList3 = Dict();
    socList4 = Dict();
    socList5 = Dict();
    for ω in Ω
      for k in fData.brList
        if fData.rateA[k] < Inf
          socList1[k,ω] = [p[k,ω],q[k,ω]];
        end
        socList2[k,ω] = [vv[k,ω],vhat[k[1],ω]/((fData.τ1[k]^2)*sqrt(2)),vhat[k[2],ω]/((fData.τ2[k]^2)*sqrt(2))];
        socList3[k,ω] = [tAux2[k,ω],tAux3[k,ω]];
        socList5[k,ω] = [wc[k,ω],ws[k,ω],vhat[k[1],ω]/((fData.τ1[k]^2)*sqrt(2)),vhat[k[2],ω]/((fData.τ2[k]^2)*sqrt(2))];
      end
      for i in fData.IDList
        socList4[i,ω] = [v[i,ω],tAux5[i,ω]];
      end
    end

    @constraint(mp,socConstraint1[k in fData.brList,ω in Ω; fData.rateA[k] < Inf],norm(socList1[k,ω]) <= fData.rateA[k]);
    @constraint(mp,socConstraint2[k in fData.brList,ω in Ω], norm(socList2[k,ω]) <= (vhat[k[1],ω]/(fData.τ1[k]^2) + vhat[k[2],ω]/(fData.τ2[k]^2))/sqrt(2));
    @constraint(mp,socConstraint3[k in fData.brList,ω in Ω],norm(socList3[k,ω]) <= tAux4[k,ω]);
    @constraint(mp,socConstraint4[i in fData.IDList,ω in Ω],norm(socList4[i,ω]) <= tAux6[i,ω]);
    @constraint(mp,socConstraint5[k in fData.brList,ω in Ω],norm(socList5[k,ω]) <= (vhat[k[1],ω]/(fData.τ1[k]^2) + vhat[k[2],ω]/(fData.τ2[k]^2))/sqrt(2));

    @constraint(mp,auxConstr2[k in fData.brList,ω in Ω],sqrt((1-cos(θu[k]))/(θu[k])^2)*((θ[k[1],ω] - fData.σ[k]) - θ[k[2],ω]) == tAux2[k,ω]);
    @constraint(mp,auxConstr3[k in fData.brList,ω in Ω],cs[k,ω] - 3/4 == tAux3[k,ω]);
    @constraint(mp,auxConstr4[k in fData.brList,ω in Ω],5/4 - cs[k,ω] == tAux4[k,ω]);
    @constraint(mp,auxConstr5[i in fData.IDList,ω in Ω],tAux5[i,ω] == vhat[i,ω] - 1/4);
    @constraint(mp,auxConstr6[i in fData.IDList,ω in Ω],tAux6[i,ω] == vhat[i,ω] + 1/4);
  else
    @constraint(mp,socConstraint1[k in fData.brList,ω in Ω; fData.rateA[k] < Inf],p[k,ω]^2 + q[k,ω]^2 <= fData.rateA[k]^2);
    @NLconstraint(mp,socConstraint2[k in fData.brList,ω in Ω],vv[k,ω]^2 <= vhat[k[1],ω]/fData.τ1[k]^2 * vhat[k[2],ω]/fData.τ2[k]^2);
    @NLconstraint(mp,socConstraint3[k in fData.brList,ω in Ω],cs[k,ω] + (1-cos(θu[k]))/(θu[k])^2*((θ[k[1],ω] - fData.σ[k]) - θ[k[2],ω])^2 <= 1);
    @NLconstraint(mp,socConstraint4[i in fData.IDList,ω in Ω],v[i,ω]^2 <= vhat[i,ω]);
    @NLconstraint(mp,socConstraint5[k in fData.brList,ω in Ω],wc[k,ω]^2 + ws[k,ω]^2 <= vhat[k[1],ω]/(fData.τ1[k]^2)*vhat[k[2],ω]/(fData.τ2[k]^2));
  end
  @constraint(mp,sinMock1[k in fData.brList,ω in Ω],
              ss[k,ω] <= cos(θu[k]/2)*((θ[k[1],ω] - fData.σ[k]) - θ[k[2],ω] - θu[k]/2) + sin(θu[k]/2));
  @constraint(mp,sinMock2[k in fData.brList,ω in Ω],
              ss[k,ω] >= cos(θu[k]/2)*((θ[k[1],ω] - fData.σ[k]) - θ[k[2],ω] + θu[k]/2) - sin(θu[k]/2));
  @constraint(mp,angleDiff1[k in fData.brList,ω in Ω],(θ[k[1],ω] - fData.σ[k]) - θ[k[2],ω] <= θDmax[k]);
  @constraint(mp,angleDiff2[k in fData.brList,ω in Ω],(θ[k[1],ω] - fData.σ[k]) - θ[k[2],ω] >= θDmin[k]);

  @constraint(mp,v2Mock2[i in fData.IDList,ω in Ω], vhat[i,ω] - (vmax[i] + vmin[i])*v[i,ω] <= -vmax[i]*vmin[i]);
  @constraint(mp,vvMock1[k in fData.brList,ω in Ω],
              vv[k,ω] >= vmin[k[1]]*v[k[2],ω]/(fData.τ1[k]*fData.τ2[k]) + vmin[k[2]]*v[k[1],ω]/(fData.τ1[k]*fData.τ2[k]) - vmin[k[1]]*vmin[k[2]]/(fData.τ1[k]*fData.τ2[k]));
  @constraint(mp,vvMock2[k in fData.brList,ω in Ω],
              vv[k,ω] >= vmax[k[1]]*v[k[2],ω]/(fData.τ1[k]*fData.τ2[k]) + vmax[k[2]]*v[k[1],ω]/(fData.τ1[k]*fData.τ2[k]) - vmax[k[1]]*vmax[k[2]]/(fData.τ1[k]*fData.τ2[k]));
  @constraint(mp,vvMock3[k in fData.brList,ω in Ω],
              vv[k,ω] <= vmin[k[1]]*v[k[2],ω]/(fData.τ1[k]*fData.τ2[k]) + vmax[k[2]]*v[k[1],ω]/(fData.τ1[k]*fData.τ2[k]) - vmin[k[1]]*vmax[k[2]]/(fData.τ1[k]*fData.τ2[k]));
  @constraint(mp,vvMock4[k in fData.brList,ω in Ω],
              vv[k,ω] <= vmax[k[1]]*v[k[2],ω]/(fData.τ1[k]*fData.τ2[k]) + vmin[k[2]]*v[k[1],ω]/(fData.τ1[k]*fData.τ2[k]) - vmax[k[1]]*vmin[k[2]]/(fData.τ1[k]*fData.τ2[k]));

  @constraint(mp,wcMock1[k in fData.brList,ω in Ω],
              wc[k,ω] >= vmin[k[1]]*vmin[k[2]]/(fData.τ1[k]*fData.τ2[k])*cs[k,ω] + csmin[k]*vv[k,ω] - vmin[k[1]]*vmin[k[2]]/(fData.τ1[k]*fData.τ2[k])*csmin[k]);
  @constraint(mp,wcMock2[k in fData.brList,ω in Ω],
              wc[k,ω] >= vmax[k[1]]*vmax[k[2]]/(fData.τ1[k]*fData.τ2[k])*cs[k,ω] + csmax[k]*vv[k,ω] - vmax[k[1]]*vmax[k[2]]/(fData.τ1[k]*fData.τ2[k])*csmax[k]);
  @constraint(mp,wcMock3[k in fData.brList,ω in Ω],
              wc[k,ω] <= vmin[k[1]]*vmin[k[2]]/(fData.τ1[k]*fData.τ2[k])*cs[k,ω] + vv[k,ω]*csmax[k] - vmin[k[1]]*vmin[k[2]]/(fData.τ1[k]*fData.τ2[k])*csmax[k]);
  @constraint(mp,wcMock4[k in fData.brList,ω in Ω],
              wc[k,ω] <= vmax[k[1]]*vmax[k[2]]/(fData.τ1[k]*fData.τ2[k])*cs[k,ω] + vv[k,ω]*csmin[k] - vmax[k[1]]*vmax[k[2]]/(fData.τ1[k]*fData.τ2[k])*csmin[k]);

  @constraint(mp,wsMock1[k in fData.brList,ω in Ω],
              ws[k,ω] >= vmin[k[1]]*vmin[k[2]]/(fData.τ1[k]*fData.τ2[k])*ss[k,ω] + ssmin[k]*vv[k,ω] - vmin[k[1]]*vmin[k[2]]/(fData.τ1[k]*fData.τ2[k])*ssmin[k]);
  @constraint(mp,wsMock2[k in fData.brList,ω in Ω],
              ws[k,ω] >= vmax[k[1]]*vmax[k[2]]/(fData.τ1[k]*fData.τ2[k])*ss[k,ω] + ssmax[k]*vv[k,ω] - vmax[k[1]]*vmax[k[2]]/(fData.τ1[k]*fData.τ2[k])*ssmax[k]);
  @constraint(mp,wsMock3[k in fData.brList,ω in Ω],
              ws[k,ω] <= vmin[k[1]]*vmin[k[2]]/(fData.τ1[k]*fData.τ2[k])*ss[k,ω] + vv[k,ω]*ssmax[k] - vmin[k[1]]*vmin[k[2]]/(fData.τ1[k]*fData.τ2[k])*ssmax[k]);
  @constraint(mp,wsMock4[k in fData.brList,ω in Ω],
              ws[k,ω] <= vmax[k[1]]*vmax[k[2]]/(fData.τ1[k]*fData.τ2[k])*ss[k,ω] + vv[k,ω]*ssmin[k] - vmax[k[1]]*vmax[k[2]]/(fData.τ1[k]*fData.τ2[k])*ssmin[k]);

  @constraint(mp,tanConstr1[k in fData.brList,ω in Ω],
                ws[k,ω] - tan(θDmax[k])*wc[k,ω] <= 0);
  @constraint(mp,tanConstr2[k in fData.brList,ω in Ω],
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

  @constraint(mp,lncConstr1[k in fData.brList,ω in Ω],
              vδ[k[1]]*vδ[k[2]]/(fData.τ1[k]*fData.τ2[k])*(wc[k,ω]*cos(θϕ[k]) + ws[k,ω]*sin(θϕ[k])) -
              vmax[k[2]]/fData.τ2[k]*cos(θδ[k])*vδ[k[2]]/fData.τ2[k]*vhat[k[1],ω]/(fData.τ1[k]^2) -
              vmax[k[1]]/fData.τ1[k]*cos(θδ[k])*vδ[k[1]]/fData.τ1[k]*vhat[k[2],ω]/(fData.τ2[k]^2) >=
              vmax[k[1]]*vmax[k[2]]/(fData.τ1[k]*fData.τ2[k])*cos(θδ[k])*(vmin[k[1]]*vmin[k[2]]/(fData.τ1[k]*fData.τ2[k]) -
              vmax[k[1]]*vmax[k[2]]/(fData.τ1[k]*fData.τ2[k])));
  @constraint(mp,lncConstr2[k in fData.brList,ω in Ω],
              vδ[k[1]]*vδ[k[2]]/(fData.τ1[k]*fData.τ2[k])*(wc[k,ω]*cos(θϕ[k]) + ws[k,ω]*sin(θϕ[k])) -
              vmin[k[2]]/fData.τ2[k]*cos(θδ[k])*vδ[k[2]]/fData.τ2[k]*vhat[k[1],ω]/(fData.τ1[k]^2) -
              vmin[k[1]]/fData.τ1[k]*cos(θδ[k])*vδ[k[1]]/fData.τ1[k]*vhat[k[2],ω]/(fData.τ2[k]^2) >=
              -vmin[k[1]]*vmin[k[2]]/(fData.τ1[k]*fData.τ2[k])*cos(θδ[k])*(vmin[k[1]]*vmin[k[2]]/(fData.τ1[k]*fData.τ2[k]) -
              vmax[k[1]]*vmax[k[2]]/(fData.τ1[k]*fData.τ2[k])));

  conWPpar1 = Dict();
  conWQpar1 = Dict();
  conWPpar2 = Dict();
  conWQpar2 = Dict();
  for i in fData.IDList
    if (abs(uData[i].WPmax) <= 1e-3)&(abs(uData[i].WPmax) > 0)
      conWPpar1[i] = 1/sqrt(abs(uData[i].WPmax));
      conWPpar2[i] = sqrt(abs(uData[i].WPmax));
    else
      conWPpar1[i] = 1;
      conWPpar2[i] = 1;
    end
    if (abs(uData[i].WQmax) <= 1e-3)&(abs(uData[i].WQmax) > 0)
      conWQpar1[i] = 1/sqrt(abs(uData[i].WQmax));
      conWQpar2[i] = sqrt(abs(uData[i].WQmax));
    else
      conWQpar1[i] = 1;
      conWQpar2[i] = 1;
    end
  end

  ypplus = Dict();
  ypminus = Dict();
  yqplus = Dict();
  yqminus = Dict();
  for ω in Ω
    for i in 1:length(fData.IDList)
      ypplus[fData.IDList[i],ω] = max(yp[ω,i]*2 - 1,0);
      ypminus[fData.IDList[i],ω] = max(1 - yp[ω,i]*2,0);
      yqplus[fData.IDList[i],ω] = max(yq[ω,i]*2 - 1,0);
      yqminus[fData.IDList[i],ω] = max(1 - yq[ω,i]*2,0);
    end
  end

  @constraint(mp,totalP[i in fData.IDList,ω in Ω], sum(p[k,ω] for k in branchDict1[i]) + vhat[i,ω]*fData.gs[i] + conWPpar2[i]*opplus[i,ω] - conWPpar2[i]*opminus[i,ω] ==
    sphatsum[i] + uData[i].GP0 + (uData[i].GPmax - uData[i].GP0)*ypplus[i,ω] + (uData[i].GPmin - uData[i].GP0)*ypminus[i,ω]);
  @constraint(mp,totalQ[i in fData.IDList,ω in Ω], sum(q[k,ω] for k in branchDict1[i]) - vhat[i,ω]*fData.bs[i] + conWQpar2[i]*oqplus[i,ω] - conWQpar2[i]*oqminus[i,ω] ==
    sqhatsum[i] + uData[i].GQ0 + (uData[i].GQmax - uData[i].GQ0)*yqplus[i,ω] + (uData[i].GQmin - uData[i].GQ0)*yqminus[i,ω]);

    @constraint(mp,opConstr1[i in fData.IDList,ω in Ω],opminus[i,ω] <= conWPpar1[i]*uData[i].WPmax);
    @constraint(mp,oqConstr1[i in fData.IDList,ω in Ω],oqminus[i,ω] <= conWQpar1[i]*uData[i].WQmax);
    @constraint(mp,opConstr2[i in fData.IDList,ω in Ω],opplus[i,ω] <= conWPpar1[i]*(uData[i].WPmax + uData[i].RESP0 +
      (uData[i].RESPmax - uData[i].RESP0)*ypplus[i,ω] + (uData[i].RESPmin - uData[i].RESP0)*ypminus[i,ω]));
    @constraint(mp,oqConstr2[i in fData.IDList,ω in Ω],oqplus[i,ω] <= conWQpar1[i]*(uData[i].WQmax + uData[i].RESQ0 +
      (uData[i].RESQmax - uData[i].RESQ0)*yqplus[i,ω] + (uData[i].RESQmin - uData[i].RESQ0)*yqminus[i,ω]));

  # append the linear cuts left
  for ω in Ω
    cLeft = [];
    for ci in 1:length(cutList)
      if !((cutList[ci][6][1][:] == yp[ω,:])&(cutList[ci][6][2][:] == yq[ω,:]))
        # remove the cut from the list
        push!(cLeft,ci)
      end
    end
    cutList = cutList[cLeft];
  end
  for ci in 1:length(cutList)
      # append the cut to mp
    @constraint(mp,cutList[ci][1] - sum(cutList[ci][2][i]*(mp[:sp][i] - cutList[ci][4][i]) for i in fData.genIDList) -
      sum(cutList[ci][3][i]*(mp[:sq][i]-cutList[ci][5][i]) for i in fData.genIDList) <= 0);
  end
  return mp;
end
