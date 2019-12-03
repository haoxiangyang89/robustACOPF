function paraSolve_One(fData,uData,yps,yqs,sphat,sqhat,vmax,vmin,θDmax,θDmin,o,numericOpt = 0,QP_opt = 0)
  # spH,spHobjo = createSH_Budget(fData,uData,yps,yqs,vmax,vmin,θDmax,θDmin);
  θu = Dict();
  for k in fData.brList
    θu[k] = max(abs(θDmax[k]),abs(θDmin[k]));
  end
  # spH = Model(solver = GurobiSolver(Threads = 1, OutputFlag = 0, NumericFocus = 3));

  if QP_opt == 0
    spH = Model(solver = GurobiSolver(Threads = 1, OutputFlag = 0, NumericFocus = numericOpt));
  else
    spH = Model(solver = IpoptSolver(print_level = 0,max_iter = 10000,linear_solver = "ma27"));
  end
  # spH = Model(solver = MosekSolver(MSK_DPAR_MIO_TOL_REL_GAP = 1e-6,MSK_DPAR_INTPNT_TOL_PFEAS = 1e-6,
  #   MSK_DPAR_MIO_TOL_FEAS = 1e-6,MSK_IPAR_LOG = 0, MSK_IPAR_NUM_THREADS = 1));
  M = 1;

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
    if abs(fData.b[k]) >= 1000
      con1Par[k] = 1/sqrt(abs(fData.b[k]));
      con2Par[k] = sqrt(abs(fData.b[k]));
    else
      con1Par[k] = 1;
      con2Par[k] = 1;
    end
  end

  # obtain the upper bound of the cs and ss terms
  csmax = Dict();
  csmin = Dict();
  ssmax = Dict();
  ssmin = Dict();
  csConst = Dict();
  ssConst = Dict();
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

  # add the strengthened bounds on ss
  for k in fData.brList
    ssmax[k] = sin(θDmax[k]);
    ssmin[k] = sin(θDmin[k]);
    csConst[k] = (cos(θDmax[k]) - cos(θDmin[k]))/(θDmax[k] - θDmin[k]);
    ssConst[k] = (sin(θDmax[k]) - sin(θDmin[k]))/(θDmax[k] - θDmin[k]);
  end

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

  # define the dual variables for the inequality A^1 x <= b^1
  # dual variables for power flow equations
  @variable(spH,λptrans[fData.brList]);
  @variable(spH,λqtrans[fData.brList]);
  @variable(spH,λps[fData.brList]);
  @variable(spH,λqs[fData.brList]);
  # dual variables for cosine relaxation
  @variable(spH,λcs1[k in fData.brList] >= 0);
  @variable(spH,λcs2[k in fData.brList] <= 0);
  @variable(spH,λcs3[k in fData.brList] <= 0);
  # dual variables for sine relaxation
  @variable(spH,λss1[k in fData.brList] >= 0);
  @variable(spH,λss2[k in fData.brList] <= 0);
  @variable(spH,λss4[k in fData.brList] >= 0);
  @variable(spH,λss5[k in fData.brList] <= 0);
  # dual variables for voltage constraints
  @variable(spH,λv[i in fData.IDList] >= 0);
  @variable(spH,λvu[i in fData.IDList] >= 0);
  @variable(spH,λvl[i in fData.IDList] <= 0);
  # dual variables for θ difference bounds
  @variable(spH,λθ1[k in fData.brList] >= 0);
  @variable(spH,λθ2[k in fData.brList] <= 0);
  @variable(spH,λθRef);
  # dual variables for vv McCormick
  @variable(spH,λvv1[k in fData.brList] <= 0);
  @variable(spH,λvv2[k in fData.brList] <= 0);
  @variable(spH,λvv3[k in fData.brList] >= 0);
  @variable(spH,λvv4[k in fData.brList] >= 0);
  # dual variables for wc McCormick
  @variable(spH,λwc1[k in fData.brList] <= 0);
  @variable(spH,λwc2[k in fData.brList] <= 0);
  @variable(spH,λwc3[k in fData.brList] >= 0);
  @variable(spH,λwc4[k in fData.brList] >= 0);
  # dual variables for ws McCormick
  @variable(spH,λws1[k in fData.brList] <= 0);
  @variable(spH,λws2[k in fData.brList] <= 0);
  @variable(spH,λws3[k in fData.brList] >= 0);
  @variable(spH,λws4[k in fData.brList] >= 0);
  # dual variables for wc/ws/cs/ss/vv/l equalities of reverse flows
  @variable(spH,λwce[k in fData.brList]);
  @variable(spH,λwse[k in fData.brList]);
  @variable(spH,λcse[k in fData.brList]);
  @variable(spH,λsse[k in fData.brList]);
  @variable(spH,λvve[k in fData.brList]);
  # dual variables for wiggle room upper bounds
  @variable(spH,λop1[i in fData.IDList] >= 0);
  @variable(spH,λoq1[i in fData.IDList] >= 0);
  @variable(spH,λop2[i in fData.IDList] >= 0);
  @variable(spH,λoq2[i in fData.IDList] >= 0);
  # dual variables for tangent constraints
  @variable(spH,λtangent1[k in fData.brList] >= 0);
  @variable(spH,λtangent2[k in fData.brList] <= 0);
  # dual variables for LNC constraints
  @variable(spH,λlnc1[k in fData.brList] <= 0);
  @variable(spH,λlnc2[k in fData.brList] <= 0);

  # define the dual variables for the inequality A^p x == s^p + u^p
  @variable(spH,-1 <= λpi[i in fData.IDList] <= 1);
  @variable(spH,-1 <= λqi[i in fData.IDList] <= 1);

  # define the dual variables for the SOC inequalities
  @variable(spH,μ1[k in fData.brList,j in 1:2; !(fData.rateA[k]==Inf)]);
  @variable(spH,μ2[k in fData.brList,j in 1:3]);
  @variable(spH,μ3[k in fData.brList,j in 1:2]);
  @variable(spH,μ4[i in fData.IDList,j in 1:2]);
  @variable(spH,μ5[k in fData.brList,j in 1:4]);
  @variable(spH,ν1[k in fData.brList; !(fData.rateA[k]==Inf)]>= 0);
  @variable(spH,ν2[k in fData.brList] >= 0);
  @variable(spH,ν3[k in fData.brList] >= 0);
  @variable(spH,ν4[i in fData.IDList]>= 0);
  @variable(spH,ν5[k in fData.brList] >= 0);

  if QP_opt == 0
    # set up the SOC constraints & linearization constraints
    for k in fData.brList
      if fData.rateA[k]<Inf
        μ1List = [μ1[k,1],μ1[k,2]];
        @constraint(spH,norm(μ1List) <= ν1[k]);
      end
      μ2List = [μ2[k,1],μ2[k,2],μ2[k,3]];
      @constraint(spH,norm(μ2List) <= ν2[k]);
      μ3List = [μ3[k,1],μ3[k,2]];
      @constraint(spH,norm(μ3List) <= ν3[k]);
      μ5List = [μ5[k,1],μ5[k,2],μ5[k,3],μ5[k,4]];
      @constraint(spH,norm(μ5List) <= ν5[k]);
    end
    for i in fData.IDList
      μ4List = [μ4[i,1],μ4[i,2]];
      @constraint(spH,norm(μ4List) <= ν4[i]);
    end
  else
    for k in fData.brList
      if fData.rateA[k]<Inf
        @constraint(spH, μ1[k,1]^2 + μ1[k,2]^2 <= ν1[k]^2);
      end
      @constraint(spH,μ2[k,1]^2+μ2[k,2]^2+μ2[k,3]^2 <= ν2[k]^2);
      @constraint(spH,μ3[k,1]^2+μ3[k,2]^2 <= ν3[k]^2);
      @constraint(spH,μ5[k,1]^2+μ5[k,2]^2+μ5[k,3]^2+μ5[k,4]^2 <= ν5[k]^2);
    end
    for i in fData.IDList
      @constraint(spH,μ4[i,1]^2 + μ4[i,2]^2 <= ν4[i]^2);
    end
  end

  # set up the x constraints
  @constraint(spH, pConstr1[k in fData.brList; fData.rateA[k]<Inf], λps[k]-μ1[k,1]+λpi[k[1]] == 0);
  @constraint(spH, pConstr2[k in fData.brList; fData.rateA[k]==Inf], λps[k]+λpi[k[1]] == 0);
  @constraint(spH, qConstr1[k in fData.brList; fData.rateA[k]<Inf], λqs[k]-μ1[k,2]+λqi[k[1]] == 0);
  @constraint(spH, qConstr2[k in fData.brList; fData.rateA[k]==Inf], λqs[k]+λqi[k[1]] == 0);
  @constraint(spH, psConstr[k in fData.brList], -λps[k]*con2Par[k] + λptrans[k] == 0);
  @constraint(spH, qsConstr[k in fData.brList], -λqs[k]*con2Par[k] + λqtrans[k] == 0);
  @constraint(spH, vConstr[i in fData.IDList], -μ4[i,1] - (vmax[i] + vmin[i])*λv[i] + λvu[i] + λvl[i] +
              sum(-vmin[k[2]]*λvv1[k]/(fData.τ1[k]*fData.τ2[k]) - vmax[k[2]]*λvv2[k]/(fData.τ1[k]*fData.τ2[k]) -
                  vmax[k[2]]*λvv3[k]/(fData.τ1[k]*fData.τ2[k]) - vmin[k[2]]*λvv4[k]/(fData.τ1[k]*fData.τ2[k]) for k in branchDict1[i]) +
              sum(-vmin[k[1]]*λvv1[k]/(fData.τ1[k]*fData.τ2[k]) - vmax[k[1]]*λvv2[k]/(fData.τ1[k]*fData.τ2[k]) -
                  vmin[k[1]]*λvv3[k]/(fData.τ1[k]*fData.τ2[k]) - vmax[k[1]]*λvv4[k]/(fData.τ1[k]*fData.τ2[k]) for k in branchDict2[i]) == 0);
  @constraint(spH,vhatConstr[i in fData.IDList], -μ4[i,2] - ν4[i] + λv[i] + λpi[i]*fData.gs[i] - λqi[i]*fData.bs[i]
              - sum(fData.g[k]*con1Par[k]*λptrans[k]/(fData.τ1[k]^2)-(fData.b[k] + fData.bc[k]/2)*con1Par[k]*λqtrans[k]/(fData.τ1[k]^2) for k in branchDict1[i])
              - sum(μ2[k,2]/(fData.τ1[k]^2*sqrt(2)) + ν2[k]/(fData.τ1[k]^2*sqrt(2)) + μ5[k,3]/(fData.τ1[k]^2*sqrt(2)) + ν5[k]/(fData.τ1[k]^2*sqrt(2)) +
              cos(θδ[k])*vδ[k[2]]/fData.τ2[k]*(vmax[k[2]]/fData.τ2[k]*λlnc1[k] + vmin[k[2]]/fData.τ2[k]*λlnc2[k])/(fData.τ1[k]^2) for k in branchDict1[i])
              - sum(μ2[k,3]/(fData.τ2[k]^2*sqrt(2)) + ν2[k]/(fData.τ2[k]^2*sqrt(2)) + μ5[k,4]/(fData.τ2[k]^2*sqrt(2)) + ν5[k]/(fData.τ2[k]^2*sqrt(2)) +
              cos(θδ[k])*vδ[k[1]]/fData.τ1[k]*(vmax[k[1]]/fData.τ1[k]*λlnc1[k] + vmin[k[1]]/fData.τ1[k]*λlnc2[k])/(fData.τ2[k]^2) for k in branchDict2[i]) == 0);
  @constraint(spH,vvConstr[k in fData.brList],
              λvv1[k]+λvv2[k]+λvv3[k]+λvv4[k] + λvve[k] - λvve[(k[2],k[1],k[3])] - μ2[k,1]
              - csmin[k]*λwc1[k] -csmax[k]*λwc2[k] - csmax[k]*λwc3[k] - csmin[k]*λwc4[k]
              - ssmin[k]*λws1[k]-ssmax[k]*λws2[k]-ssmax[k]*λws3[k]-ssmin[k]*λws4[k] == 0);
  @constraint(spH,wcConstr[k in fData.brList],
              fData.g[k]*con1Par[k]*λptrans[k] - fData.b[k]*con1Par[k]*λqtrans[k] - μ5[k,1] +
              λwc1[k] + λwc2[k] + λwc3[k] + λwc4[k] + λwce[k] - λwce[(k[2],k[1],k[3])] -
              λtangent1[k]*tan(θDmax[k]) - λtangent2[k]*tan(θDmin[k]) +
              vδ[k[1]]*vδ[k[2]]/(fData.τ1[k]*fData.τ2[k])*cos(θϕ[k])*(λlnc1[k] + λlnc2[k]) == 0);
  @constraint(spH,wsConstr[k in fData.brList],
              fData.b[k]*con1Par[k]*λptrans[k]+fData.g[k]*con1Par[k]*λqtrans[k] + λtangent1[k] + λtangent2[k] - μ5[k,2] +
              λws1[k] + λws2[k] + λws3[k] + λws4[k] + λwse[k] + λwse[(k[2],k[1],k[3])] +
              vδ[k[1]]*vδ[k[2]]/(fData.τ1[k]*fData.τ2[k])*sin(θϕ[k])*(λlnc1[k] + λlnc2[k]) == 0);
  @constraint(spH,csConstr[k in fData.brList],
              λcs1[k] + λcs2[k] + λcs3[k] - μ3[k,2] + ν3[k] + λcse[k] - λcse[(k[2],k[1],k[3])]
              - vmin[k[1]]*vmin[k[2]]/(fData.τ1[k]*fData.τ2[k])*λwc1[k] - vmax[k[1]]*vmax[k[2]]*λwc2[k]/(fData.τ1[k]*fData.τ2[k])
              - vmin[k[1]]*vmin[k[2]]/(fData.τ1[k]*fData.τ2[k])*λwc3[k] - vmax[k[1]]*vmax[k[2]]*λwc4[k]/(fData.τ1[k]*fData.τ2[k]) == 0);
  @constraint(spH,ssConstr[k in fData.brList],
              λss1[k] + λss2[k] + λss4[k] + λss5[k] + λsse[k] + λsse[(k[2],k[1],k[3])]
              - vmin[k[1]]*vmin[k[2]]*λws1[k]/(fData.τ1[k]*fData.τ2[k]) - vmax[k[1]]*vmax[k[2]]*λws2[k]/(fData.τ1[k]*fData.τ2[k])
              - vmin[k[1]]*vmin[k[2]]*λws3[k]/(fData.τ1[k]*fData.τ2[k]) - vmax[k[1]]*vmax[k[2]]*λws4[k]/(fData.τ1[k]*fData.τ2[k]) == 0);
  refBusInd = Dict();
  for i in fData.IDList
    if fData.bType[i] == 3
      refBusInd[i] = 1;
    else
      refBusInd[i] = 0;
    end
  end
  @constraint(spH,θConstr[i in fData.IDList],
              λθRef*refBusInd[i] + sum(λθ1[k] + λθ2[k] - csConst[k]*λcs3[k] - cos(θu[k]/2)*λss4[k] - cos(θu[k]/2)*λss5[k]
              - sqrt(1-cos(θu[k]))/θu[k]*μ3[k,1] for k in branchDict1[i])
              + sum(-λθ1[k] - λθ2[k] + csConst[k]*λcs3[k] + cos(θu[k]/2)*λss4[k] + cos(θu[k]/2)*λss5[k]
              + sqrt(1-cos(θu[k]))/(θu[k])*μ3[k,1] for k in branchDict2[i]) == 0);
  @constraint(spH,opConstr1[i in fData.IDList],λop1[i] >= λpi[i]);
  @constraint(spH,opConstr2[i in fData.IDList],λop2[i] >= -λpi[i]);
  @constraint(spH,oqConstr1[i in fData.IDList],λoq1[i] >= λqi[i]);
  @constraint(spH,oqConstr2[i in fData.IDList],λoq2[i] >= -λqi[i]);

  # set up the objective function
  ypplus = Dict();
  ypminus = Dict();
  yqplus = Dict();
  yqminus = Dict();
  for i in 1:length(fData.IDList)
    ypplus[fData.IDList[i]] = max(yps[i]*2 - 1,0);
    ypminus[fData.IDList[i]] = max(1 - yps[i]*2,0);
    yqplus[fData.IDList[i]] = max(yqs[i]*2 - 1,0);
    yqminus[fData.IDList[i]] = max(1 - yqs[i]*2,0);
  end

  # set up the objective function
  @expression(spH,restTerm,-sum(fData.rateA[k]*ν1[k] for k in fData.brList if !(fData.rateA[k] == Inf)) -
            sum(ypplus[i]*λpi[i]*(uData[i].GPmax - uData[i].GP0) + ypminus[i]*λpi[i]*(-uData[i].GP0 + uData[i].GPmin) + uData[i].GP0*λpi[i] +
            yqplus[i]*λqi[i]*(uData[i].GQmax - uData[i].GQ0) + yqminus[i]*λqi[i]*(-uData[i].GQ0 + uData[i].GQmin) + uData[i].GQ0*λqi[i] +
            -μ4[i,2]/4 + ν4[i]/4 + λvu[i]*vmax[i] + λvl[i]*vmin[i] + λop1[i]*uData[i].WPmax + λoq1[i]*uData[i].WQmax +
            λop2[i]*(uData[i].WPmax + uData[i].RESP0 + ypplus[i]*(uData[i].RESPmax - uData[i].RESP0) + ypminus[i]*(uData[i].RESPmin - uData[i].RESP0)) +
            λoq2[i]*(uData[i].WQmax + uData[i].RESQ0 + yqplus[i]*(uData[i].RESQmax - uData[i].RESQ0) + yqminus[i]*(uData[i].RESQmin - uData[i].RESQ0)) - vmax[i]*vmin[i]*λv[i] for i in fData.IDList) -
            sum(-3/4*μ3[k,2] - sqrt((1-cos(θu[k]))/(θu[k])^2)*fData.σ[k]*μ3[k,1] + 5/4*ν3[k] + csmax[k]*λcs1[k] + csmin[k]*λcs2[k] + (cos(θDmin[k]) - csConst[k]*(fData.σ[k] + θDmin[k]))*λcs3[k]
              + (sin(θu[k]/2) - θu[k]/2*cos(θu[k]/2))*(λss4[k] - λss5[k]) + ssmax[k]*λss1[k] + ssmin[k]*λss2[k]
              + (θDmax[k] + fData.σ[k])*λθ1[k] + (θDmin[k] + fData.σ[k])*λθ2[k] - fData.σ[k]*(λss4[k] + λss5[k])
              + vmax[k[1]]*vmax[k[2]]/(fData.τ1[k]*fData.τ2[k])*cos(θδ[k])*(vmin[k[1]]*vmin[k[2]] - vmax[k[1]]*vmax[k[2]])/(fData.τ1[k]*fData.τ2[k])*λlnc1[k]
              - vmin[k[1]]*vmin[k[2]]/(fData.τ1[k]*fData.τ2[k])*cos(θδ[k])*(vmin[k[1]]*vmin[k[2]] - vmax[k[1]]*vmax[k[2]])/(fData.τ1[k]*fData.τ2[k])*λlnc2[k]
              - vmin[k[1]]*vmin[k[2]]/(fData.τ1[k]*fData.τ2[k])*λvv1[k] - vmax[k[1]]*vmax[k[2]]/(fData.τ1[k]*fData.τ2[k])*λvv2[k]
              - vmin[k[1]]*vmax[k[2]]/(fData.τ1[k]*fData.τ2[k])*λvv3[k] - vmax[k[1]]*vmin[k[2]]/(fData.τ1[k]*fData.τ2[k])*λvv4[k]
              - csmin[k]*vmin[k[1]]*vmin[k[2]]/(fData.τ1[k]*fData.τ2[k])*λwc1[k] - csmax[k]*vmax[k[1]]*vmax[k[2]]/(fData.τ1[k]*fData.τ2[k])*λwc2[k]
              - csmax[k]*vmin[k[1]]*vmin[k[2]]/(fData.τ1[k]*fData.τ2[k])*λwc3[k] - csmin[k]*vmax[k[1]]*vmax[k[2]]/(fData.τ1[k]*fData.τ2[k])*λwc4[k]
              - ssmin[k]*vmin[k[1]]*vmin[k[2]]/(fData.τ1[k]*fData.τ2[k])*λws1[k] - ssmax[k]*vmax[k[1]]*vmax[k[2]]/(fData.τ1[k]*fData.τ2[k])*λws2[k]
              - ssmax[k]*vmin[k[1]]*vmin[k[2]]/(fData.τ1[k]*fData.τ2[k])*λws3[k] - ssmin[k]*vmax[k[1]]*vmax[k[2]]/(fData.τ1[k]*fData.τ2[k])*λws4[k]
              for k in fData.brList));
  for i in fData.genIDList
    # update the subproblem objective function with the current master solutions
    restTerm -= (sphat[i]*spH[:λpi][fData.Loc[i]] + sqhat[i]*spH[:λqi][fData.Loc[i]]);
  end
  @objective(spH,Max,restTerm);
  solve(spH);
  λpr = zeros(length(fData.genIDList));
  λqr = zeros(length(fData.genIDList));
  for i in fData.genIDList
    λpr[i] = getvalue(spH[:λpi][fData.Loc[i]]);
    λqr[i] = getvalue(spH[:λqi][fData.Loc[i]]);
  end
  return λpr,λqr,getobjectivevalue(spH),o;
end

function paraSolve_Multiple(fData,uData,Ωs,yps,yqs,sphat,sqhat,vmax,vmin,θDmax,θDmin,numericOpt = 0)
  spInfoList = [];
  spInfoList = pmap(ω -> paraSolve_One(fData,uData,yps[ω,:],yqs[ω,:],sphat,sqhat,vmax,vmin,θDmax,θDmin,ω,numericOpt), Ωs);
  return spInfoList;
end
