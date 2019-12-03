# function to generate the voltage bounds, including transformers
function solveBV_SF(fData,uData,ibus,mode,vmax,vmin,θDmax,θDmin,QP_opt = 1)
  θu = Dict();
  for k in fData.brList
    θu[k] = max(abs(θDmax[k]),abs(θDmin[k]));
  end
  # scaling issue exists for this primal problem.
  if QP_opt == 1
    mp = Model(solver = IpoptSolver(print_level = 0,max_iter = 10000,linear_solver = "ma27"));
  else
    mp = Model(solver = GurobiSolver(NumericFocus = 3, OutputFlag = 0, Threads = 1));
  end
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

  @variable(mp, fData.Pmin[i] <= sp[i in fData.genIDList] <= fData.Pmax[i]);
  @variable(mp, fData.Qmin[i] <= sq[i in fData.genIDList] <= fData.Qmax[i]);

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
  # set up the reference bus
  for i in fData.IDList
    if fData.bType[i] == 3
      @constraint(mp,θ[i] == 0);
    end
  end

  @variable(mp,tAux1[k in fData.brList] >= 0);
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
  # @constraint(mp,lEquality[k in fData.brList;k[1] < k[2]],l[k] == l[(k[2],k[1],k[3])]);
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

  if QP_opt == 1
    @constraint(mp,socConstraint1[k in fData.brList; fData.rateA[k] < Inf],p[k]^2 + q[k]^2 <= fData.rateA[k]^2);
    @NLconstraint(mp,socConstraint2[k in fData.brList],vv[k]^2 <= vhat[k[1]]/fData.τ1[k]^2 * vhat[k[2]]/fData.τ2[k]^2);
    @NLconstraint(mp,socConstraint3[k in fData.brList],cs[k] + (1-cos(θu[k]))/(θu[k])^2*((θ[k[1]] - fData.σ[k]) - θ[k[2]])^2 <= 1);
    @NLconstraint(mp,socConstraint4[i in fData.IDList],v[i]^2 <= vhat[i]);
    @NLconstraint(mp,socConstraint5[k in fData.brList],wc[k]^2 + ws[k]^2 <= vhat[k[1]]/(fData.τ1[k]^2)*vhat[k[2]]/(fData.τ2[k]^2));
  else
    socList1 = Dict();
    socList3 = Dict();
    socList4 = Dict();
    socList5 = Dict();
    for k in fData.brList
      if fData.rateA[k] < Inf
        socList1[k] = [p[k],q[k]];
      end
      socList3[k] = [tAux2[k],tAux3[k]];
      socList5[k] = [wc[k],ws[k],vhat[k[1]]/((fData.τ1[k]^2)*sqrt(2)),vhat[k[2]]/((fData.τ2[k]^2)*sqrt(2))];
    end
    for i in fData.IDList
      socList4[i] = [v[i],tAux5[i]];
    end

    @constraint(mp,socConstraint1[k in fData.brList; fData.rateA[k] < Inf],norm(socList1[k]) <= fData.rateA[k]);
    @constraint(mp,socConstraint3[k in fData.brList],norm(socList3[k]) <= tAux4[k]);
    @constraint(mp,socConstraint4[i in fData.IDList],norm(socList4[i]) <= tAux6[i]);
    @constraint(mp,socConstraint5[k in fData.brList],norm(socList5[k]) <= (vhat[k[1]]/(fData.τ1[k]^2) + vhat[k[2]]/(fData.τ2[k]^2))/sqrt(2));

    @constraint(mp,auxConstr2[k in fData.brList],sqrt((1-cos(θu[k]))/(θu[k])^2)*((θ[k[1]] - fData.σ[k]) - θ[k[2]]) == tAux2[k]);
    @constraint(mp,auxConstr3[k in fData.brList],cs[k] - 3/4 == tAux3[k]);
    @constraint(mp,auxConstr4[k in fData.brList],5/4 - cs[k] == tAux4[k]);
    @constraint(mp,auxConstr5[i in fData.IDList],tAux5[i] == vhat[i] - 1/4);
    @constraint(mp,auxConstr6[i in fData.IDList],tAux6[i] == vhat[i] + 1/4);
  end
  @constraint(mp,sinMock1[k in fData.brList],
              ss[k] <= cos(θu[k]/2)*((θ[k[1]] - fData.σ[k]) - θ[k[2]] - θu[k]/2) + sin(θu[k]/2));
  @constraint(mp,sinMock2[k in fData.brList],
              ss[k] >= cos(θu[k]/2)*((θ[k[1]] - fData.σ[k]) - θ[k[2]] + θu[k]/2) - sin(θu[k]/2));
  @constraint(mp,v2Mock2[i in fData.IDList], vhat[i] - (vmax[i] + vmin[i])*v[i] <= -vmax[i]*vmin[i]);
  @constraint(mp,angleDiff1[k in fData.brList],(θ[k[1]] - fData.σ[k]) - θ[k[2]] <= θDmax[k]);
  @constraint(mp,angleDiff2[k in fData.brList],(θ[k[1]] - fData.σ[k]) - θ[k[2]] >= θDmin[k]);
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

  sphatsum = Dict();
  for i in fData.IDList
    sphatsum[i] = @expression(mp,0);
    if i in keys(fData.LocRev)
      for j in fData.LocRev[i]
        sphatsum[i] += sp[j];
      end
    end
  end

  sqhatsum = Dict();
  for i in fData.IDList
    sqhatsum[i] = @expression(mp,0);
    if i in keys(fData.LocRev)
      for j in fData.LocRev[i]
        sqhatsum[i] += sq[j];
      end
    end
  end

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

  @constraint(mp,opConstr1[i in fData.IDList],opminus[i] <= conWQpar1[i]*uData[i].WPmax);
  @constraint(mp,oqConstr1[i in fData.IDList],oqminus[i] <= conWQpar1[i]*uData[i].WQmax);
  @constraint(mp,opConstr2[i in fData.IDList],opplus[i] <= conWQpar1[i]*(uData[i].WPmax + uData[i].RESPmax));
  @constraint(mp,oqConstr2[i in fData.IDList],oqplus[i] <= conWQpar1[i]*(uData[i].WQmax + uData[i].RESQmax));

  @variable(mp,uData[i].GPmin - 1e-7 <= up[i in fData.IDList] <= uData[i].GPmax + 1e-7);
  @variable(mp,uData[i].GQmin - 1e-7 <= uq[i in fData.IDList] <= uData[i].GQmax + 1e-7);

  @constraint(mp,totalP[i in fData.IDList], sum(p[k] for k in branchDict1[i]) + vhat[i]*fData.gs[i]
    + opplus[i]*conWPpar2[i] - opminus[i]*conWPpar2[i] == sphatsum[i] + up[i]);
  @constraint(mp,totalQ[i in fData.IDList], sum(q[k] for k in branchDict1[i]) - vhat[i]*fData.bs[i]
    + oqplus[i]*conWQpar2[i] - oqminus[i]*conWQpar2[i] == sqhatsum[i] + uq[i]);

  if mode == 1
    @objective(mp, Max, v[ibus]);
  else
    @objective(mp, Min, v[ibus]);
  end
  mpstatus = solve(mp);
  if mpstatus == :Optimal
    optmp = getobjectivevalue(mp);
  else
    if mode == 1
      optmp = vmax[ibus];
    else
      optmp = vmin[ibus];
    end
  end
  return optmp;
end

# function to generate the phaseangle bounds
function solveBtheta_SF(fData,uData,ibranch,mode,vmax,vmin,θDmax,θDmin,QP_opt = 1)
  θu = Dict();
  for k in fData.brList
    θu[k] = max(abs(θDmax[k]),abs(θDmin[k]));
  end
  # scaling issue exists for this primal problem.
  if QP_opt == 1
    mp = Model(solver = IpoptSolver(print_level = 0,max_iter = 10000,linear_solver = "ma27"))
  else
    mp = Model(solver = GurobiSolver(NumericFocus = 3, OutputFlag = 0, Threads = 1));
  end
  #mp = Model(solver = MosekSolver(MSK_DPAR_MIO_TOL_REL_GAP = 1e-7,MSK_DPAR_INTPNT_TOL_PFEAS = 1e-9,MSK_DPAR_MIO_TOL_FEAS = 1e-9,MSK_IPAR_LOG = 0));

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

  @variable(mp, fData.Pmin[i] <= sp[i in fData.genIDList] <= fData.Pmax[i]);
  @variable(mp, fData.Qmin[i] <= sq[i in fData.genIDList] <= fData.Qmax[i]);

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
  # set up the reference bus
  for i in fData.IDList
    if fData.bType[i] == 3
      @constraint(mp,θ[i] == 0);
    end
  end

  @variable(mp,tAux1[k in fData.brList] >= 0);
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
  # @constraint(mp,lEquality[k in fData.brList;k[1] < k[2]],l[k] == l[(k[2],k[1],k[3])]);
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

  if QP_opt == 1
    @constraint(mp,socConstraint1[k in fData.brList; fData.rateA[k] < Inf],p[k]^2 + q[k]^2 <= fData.rateA[k]^2);
    @NLconstraint(mp,socConstraint2[k in fData.brList],vv[k]^2 <= vhat[k[1]]/fData.τ1[k]^2 * vhat[k[2]]/fData.τ2[k]^2);
    @NLconstraint(mp,socConstraint3[k in fData.brList],cs[k] + (1-cos(θu[k]))/(θu[k])^2*((θ[k[1]] - fData.σ[k]) - θ[k[2]])^2 <= 1);
    @NLconstraint(mp,socConstraint4[i in fData.IDList],v[i]^2 <= vhat[i]);
    @NLconstraint(mp,socConstraint5[k in fData.brList],wc[k]^2 + ws[k]^2 <= vhat[k[1]]/(fData.τ1[k]^2)*vhat[k[2]]/(fData.τ2[k]^2));
  else
    socList1 = Dict();
    # socList2 = Dict();
    socList3 = Dict();
    socList4 = Dict();
    socList5 = Dict();
    for k in fData.brList
      if fData.rateA[k] < Inf
        socList1[k] = [p[k],q[k]];
      end
      socList3[k] = [tAux2[k],tAux3[k]];
      socList5[k] = [wc[k],ws[k],vhat[k[1]]/((fData.τ1[k]^2)*sqrt(2)),vhat[k[2]]/((fData.τ2[k]^2)*sqrt(2))];
    end
    for i in fData.IDList
      socList4[i] = [v[i],tAux5[i]];
    end

    @constraint(mp,socConstraint1[k in fData.brList; fData.rateA[k] < Inf],norm(socList1[k]) <= fData.rateA[k]);
    @constraint(mp,socConstraint3[k in fData.brList],norm(socList3[k]) <= tAux4[k]);
    @constraint(mp,socConstraint4[i in fData.IDList],norm(socList4[i]) <= tAux6[i]);
    @constraint(mp,socConstraint5[k in fData.brList],norm(socList5[k]) <= (vhat[k[1]]/(fData.τ1[k]^2) + vhat[k[2]]/(fData.τ2[k]^2))/sqrt(2));

    @constraint(mp,auxConstr2[k in fData.brList],sqrt((1-cos(θu[k]))/(θu[k])^2)*((θ[k[1]] - fData.σ[k]) - θ[k[2]]) == tAux2[k]);
    @constraint(mp,auxConstr3[k in fData.brList],cs[k] - 3/4 == tAux3[k]);
    @constraint(mp,auxConstr4[k in fData.brList],5/4 - cs[k] == tAux4[k]);
    @constraint(mp,auxConstr5[i in fData.IDList],tAux5[i] == vhat[i] - 1/4);
    @constraint(mp,auxConstr6[i in fData.IDList],tAux6[i] == vhat[i] + 1/4);
  end

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

  sphatsum = Dict();
  for i in fData.IDList
    sphatsum[i] = @expression(mp,0);
    if i in keys(fData.LocRev)
      for j in fData.LocRev[i]
        sphatsum[i] += sp[j];
      end
    end
  end

  sqhatsum = Dict();
  for i in fData.IDList
    sqhatsum[i] = @expression(mp,0);
    if i in keys(fData.LocRev)
      for j in fData.LocRev[i]
        sqhatsum[i] += sq[j];
      end
    end
  end

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

  @constraint(mp,opConstr1[i in fData.IDList],opminus[i] <= conWQpar1[i]*uData[i].WPmax);
  @constraint(mp,oqConstr1[i in fData.IDList],oqminus[i] <= conWQpar1[i]*uData[i].WQmax);
  @constraint(mp,opConstr2[i in fData.IDList],opplus[i] <= conWQpar1[i]*(uData[i].WPmax + uData[i].RESPmax));
  @constraint(mp,oqConstr2[i in fData.IDList],oqplus[i] <= conWQpar1[i]*(uData[i].WQmax + uData[i].RESQmax));

  @variable(mp,uData[i].GPmin - 1e-7 <= up[i in fData.IDList] <= uData[i].GPmax + 1e-7);
  @variable(mp,uData[i].GQmin - 1e-7 <= uq[i in fData.IDList] <= uData[i].GQmax + 1e-7);

  @constraint(mp,totalP[i in fData.IDList], sum(p[k] for k in branchDict1[i]) + vhat[i]*fData.gs[i]
    + opplus[i]*conWPpar2[i] - opminus[i]*conWPpar2[i] == sphatsum[i] + up[i]);
  @constraint(mp,totalQ[i in fData.IDList], sum(q[k] for k in branchDict1[i]) - vhat[i]*fData.bs[i]
    + oqplus[i]*conWQpar2[i] - oqminus[i]*conWQpar2[i] == sqhatsum[i] + uq[i]);

  if mode == 1
    @objective(mp, Max, θ[ibranch[1]] - θ[ibranch[2]]);
  else
    @objective(mp, Min, θ[ibranch[1]] - θ[ibranch[2]]);
  end
  mpstatus = solve(mp);
  if mpstatus == :Optimal
    optmp = getobjectivevalue(mp);
  else
    if mode == 1
      optmp = θDmax[ibranch] + fData.σ[ibranch];
    else
      optmp = θDmin[ibranch] + fData.σ[ibranch];
    end
  end
  return optmp - fData.σ[ibranch];
end
