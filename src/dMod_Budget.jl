# allow wiggle room in the second stage, including transformers

# this is the function to create the initial master problem
# remained unchanged from initMod.jl

# we allow some wiggle room in the second stage, which will change the subproblem
function createS_Budget(fData::fixedData,uData::Dict{Any,Any},vmax,vmin,θDmax,θDmin,Γ = 1,fixedy = [],numericOpt = 0,PresolveOpt = -1)
  # initialize the sub problem
  # obtian θu from the θDmax and θDmin
  θu = Dict();
  for k in fData.brList
    θu[k] = max(abs(θDmax[k]),abs(θDmin[k]));
  end
  # input: network parameters from the data variable
  sprob = Model(solver = GurobiSolver(Presolve = PresolveOpt,NumericFocus = numericOpt, OutputFlag = 0, Threads = 20));
  # sprob = Model(solver = MosekSolver(MSK_DPAR_MIO_TOL_REL_GAP = 1e-6,MSK_DPAR_INTPNT_TOL_PFEAS = 1e-6,
  #   MSK_DPAR_MIO_TOL_FEAS = 1e-6,MSK_IPAR_LOG = 0, MSK_IPAR_NUM_THREADS = 20));
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
  @variable(sprob,λptrans[fData.brList]);
  @variable(sprob,λqtrans[fData.brList]);
  @variable(sprob,λps[fData.brList]);
  @variable(sprob,λqs[fData.brList]);
  # dual variables for cosine relaxation
  @variable(sprob,λcs1[k in fData.brList] >= 0);
  @variable(sprob,λcs2[k in fData.brList] <= 0);
  @variable(sprob,λcs3[k in fData.brList] <= 0);
  # dual variables for sine relaxation
  @variable(sprob,λss1[k in fData.brList] >= 0);
  @variable(sprob,λss2[k in fData.brList] <= 0);
  @variable(sprob,λss4[k in fData.brList] >= 0);
  @variable(sprob,λss5[k in fData.brList] <= 0);
  # dual variables for voltage constraints
  @variable(sprob,λv[i in fData.IDList] >= 0);
  @variable(sprob,λvu[i in fData.IDList] >= 0);
  @variable(sprob,λvl[i in fData.IDList] <= 0);
  # dual variables for θ difference bounds
  @variable(sprob,λθ1[k in fData.brList] >= 0);
  @variable(sprob,λθ2[k in fData.brList] <= 0);
  @variable(sprob,λθRef);
  # dual variables for vv McCormick
  @variable(sprob,λvv1[k in fData.brList] <= 0);
  @variable(sprob,λvv2[k in fData.brList] <= 0);
  @variable(sprob,λvv3[k in fData.brList] >= 0);
  @variable(sprob,λvv4[k in fData.brList] >= 0);
  # dual variables for wc McCormick
  @variable(sprob,λwc1[k in fData.brList] <= 0);
  @variable(sprob,λwc2[k in fData.brList] <= 0);
  @variable(sprob,λwc3[k in fData.brList] >= 0);
  @variable(sprob,λwc4[k in fData.brList] >= 0);
  # dual variables for ws McCormick
  @variable(sprob,λws1[k in fData.brList] <= 0);
  @variable(sprob,λws2[k in fData.brList] <= 0);
  @variable(sprob,λws3[k in fData.brList] >= 0);
  @variable(sprob,λws4[k in fData.brList] >= 0);
  # dual variables for wc/ws/cs/ss/vv/l equalities of reverse flows
  @variable(sprob,λwce[k in fData.brList]);
  @variable(sprob,λwse[k in fData.brList]);
  @variable(sprob,λcse[k in fData.brList]);
  @variable(sprob,λsse[k in fData.brList]);
  @variable(sprob,λvve[k in fData.brList]);
  # dual variables for wiggle room upper bounds
  @variable(sprob,λop1[i in fData.IDList] >= 0);
  @variable(sprob,λop2[i in fData.IDList] >= 0);
  @variable(sprob,λoq1[i in fData.IDList] >= 0);
  @variable(sprob,λoq2[i in fData.IDList] >= 0);
  # dual variables for tangent constraints
  @variable(sprob,λtangent1[k in fData.brList] >= 0);
  @variable(sprob,λtangent2[k in fData.brList] <= 0);
  # dual variables for LNC constraints
  @variable(sprob,λlnc1[k in fData.brList] <= 0);
  @variable(sprob,λlnc2[k in fData.brList] <= 0);

  # define the dual variables for the inequality A^p x == s^p + u^p
  @variable(sprob,-1 <= λpi[i in fData.IDList] <= 1);
  @variable(sprob,-1 <= λqi[i in fData.IDList] <= 1);

  # define the dual variables for the SOC inequalities
  @variable(sprob,μ1[k in fData.brList,j in 1:2; !(fData.rateA[k]==Inf)]);
  @variable(sprob,μ2[k in fData.brList,j in 1:3]);
  @variable(sprob,μ3[k in fData.brList,j in 1:2]);
  @variable(sprob,μ4[i in fData.IDList,j in 1:2]);
  @variable(sprob,μ5[k in fData.brList,j in 1:4]);
  @variable(sprob,ν1[k in fData.brList; !(fData.rateA[k]==Inf)]>= 0);
  @variable(sprob,ν2[k in fData.brList] >= 0);
  @variable(sprob,ν3[k in fData.brList] >= 0);
  @variable(sprob,ν4[i in fData.IDList]>= 0);
  @variable(sprob,ν5[k in fData.brList] >= 0);

  # define the auxiliary variables: y & r

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
  @variable(sprob, ypplus[i in fData.IDList], Bin);
  @variable(sprob, rpplus[i in fData.IDList;!(i in rmListP)]);
  @variable(sprob, ypminus[i in fData.IDList], Bin);
  @variable(sprob, rpminus[i in fData.IDList;!(i in rmListP)]);
  @variable(sprob, yqplus[i in fData.IDList], Bin);
  @variable(sprob, rqplus[i in fData.IDList;!(i in rmListQ)]);
  @variable(sprob, yqminus[i in fData.IDList], Bin);
  @variable(sprob, rqminus[i in fData.IDList;!(i in rmListQ)]);

  # set up the SOC constraints & linearization constraints
  for k in fData.brList
    if fData.rateA[k]<Inf
      μ1List = [μ1[k,1],μ1[k,2]];
      @constraint(sprob,norm(μ1List) <= ν1[k]);
    end
    μ2List = [μ2[k,1],μ2[k,2],μ2[k,3]];
    @constraint(sprob,norm(μ2List) <= ν2[k]);
    μ3List = [μ3[k,1],μ3[k,2]];
    @constraint(sprob,norm(μ3List) <= ν3[k]);
    μ5List = [μ5[k,1],μ5[k,2],μ5[k,3],μ5[k,4]];
    @constraint(sprob,norm(μ5List) <= ν5[k]);
  end
  for i in fData.IDList
    μ4List = [μ4[i,1],μ4[i,2]];
    @constraint(sprob,norm(μ4List) <= ν4[i]);
  end
  @constraint(sprob, ubRpplus1[i in fData.IDList;!(i in rmListP)], rpplus[i] <= ypplus[i]);
  @constraint(sprob, lbRpplus1[i in fData.IDList;!(i in rmListP)], rpplus[i] >= -ypplus[i]);
  @constraint(sprob, ubRpplus2[i in fData.IDList;!(i in rmListP)], rpplus[i] <= λpi[i] + (1-ypplus[i])*M);
  @constraint(sprob, lbRpplus2[i in fData.IDList;!(i in rmListP)], rpplus[i] >= λpi[i] - (1-ypplus[i])*M);
  @constraint(sprob, ubRpminus1[i in fData.IDList;!(i in rmListP)], rpminus[i] <= ypminus[i]);
  @constraint(sprob, lbRpminus1[i in fData.IDList;!(i in rmListP)], rpminus[i] >= -ypminus[i]);
  @constraint(sprob, ubRpminus2[i in fData.IDList;!(i in rmListP)], rpminus[i] <= λpi[i] + (1-ypminus[i])*M);
  @constraint(sprob, lbRpminus2[i in fData.IDList;!(i in rmListP)], rpminus[i] >= λpi[i] - (1-ypminus[i])*M);
  @constraint(sprob, ubRqplus1[i in fData.IDList;!(i in rmListQ)], rqplus[i] <= yqplus[i]);
  @constraint(sprob, lbRqplus1[i in fData.IDList;!(i in rmListQ)], rqplus[i] >= -yqplus[i]);
  @constraint(sprob, ubRqplus2[i in fData.IDList;!(i in rmListQ)], rqplus[i] <= λqi[i] + (1-yqplus[i])*M);
  @constraint(sprob, lbRqplus2[i in fData.IDList;!(i in rmListQ)], rqplus[i] >= λqi[i] - (1-yqplus[i])*M);
  @constraint(sprob, ubRqminus1[i in fData.IDList;!(i in rmListQ)], rqminus[i] <= yqminus[i]);
  @constraint(sprob, lbRqminus1[i in fData.IDList;!(i in rmListQ)], rqminus[i] >= -yqminus[i]);
  @constraint(sprob, ubRqminus2[i in fData.IDList;!(i in rmListQ)], rqminus[i] <= λqi[i] + (1-yqminus[i])*M);
  @constraint(sprob, lbRqminus2[i in fData.IDList;!(i in rmListQ)], rqminus[i] >= λqi[i] - (1-yqminus[i])*M);
  @constraint(sprob, oneVertexp[i in fData.IDList;!(i in rmListP)], ypplus[i]+ypminus[i] <= 1);
  @constraint(sprob, oneVertexq[i in fData.IDList;!(i in rmListQ)], yqplus[i]+yqminus[i] <= 1);
  @constraint(sprob, budget, sum(ypplus[i]+ypminus[i] for i in fData.IDList if !(i in rmListP)) +
    sum(yqplus[i] + yqminus[i] for i in fData.IDList if !(i in rmListQ)) <= Γ);

  # set up the x constraints
  @constraint(sprob, pConstr1[k in fData.brList; fData.rateA[k]<Inf], λps[k]-μ1[k,1]+λpi[k[1]] == 0);
  @constraint(sprob, pConstr2[k in fData.brList; fData.rateA[k]==Inf], λps[k]+λpi[k[1]] == 0);
  @constraint(sprob, qConstr1[k in fData.brList; fData.rateA[k]<Inf], λqs[k]-μ1[k,2]+λqi[k[1]] == 0);
  @constraint(sprob, qConstr2[k in fData.brList; fData.rateA[k]==Inf], λqs[k]+λqi[k[1]] == 0);
  @constraint(sprob, psConstr[k in fData.brList], -λps[k]*con2Par[k] + λptrans[k] == 0);
  @constraint(sprob, qsConstr[k in fData.brList], -λqs[k]*con2Par[k] + λqtrans[k] == 0);
  @constraint(sprob, vConstr[i in fData.IDList], -μ4[i,1] - (vmax[i] + vmin[i])*λv[i] + λvu[i] + λvl[i] +
              sum(-vmin[k[2]]*λvv1[k]/(fData.τ1[k]*fData.τ2[k]) - vmax[k[2]]*λvv2[k]/(fData.τ1[k]*fData.τ2[k]) -
                  vmax[k[2]]*λvv3[k]/(fData.τ1[k]*fData.τ2[k]) - vmin[k[2]]*λvv4[k]/(fData.τ1[k]*fData.τ2[k]) for k in branchDict1[i]) +
              sum(-vmin[k[1]]*λvv1[k]/(fData.τ1[k]*fData.τ2[k]) - vmax[k[1]]*λvv2[k]/(fData.τ1[k]*fData.τ2[k]) -
                  vmin[k[1]]*λvv3[k]/(fData.τ1[k]*fData.τ2[k]) - vmax[k[1]]*λvv4[k]/(fData.τ1[k]*fData.τ2[k]) for k in branchDict2[i]) == 0);
  @constraint(sprob,vhatConstr[i in fData.IDList], -μ4[i,2] - ν4[i] + λv[i] + λpi[i]*fData.gs[i] - λqi[i]*fData.bs[i]
              - sum(fData.g[k]*con1Par[k]*λptrans[k]/(fData.τ1[k]^2)-(fData.b[k] + fData.bc[k]/2)*con1Par[k]*λqtrans[k]/(fData.τ1[k]^2) for k in branchDict1[i])
              - sum(μ2[k,2]/(fData.τ1[k]^2*sqrt(2)) + ν2[k]/(fData.τ1[k]^2*sqrt(2)) + μ5[k,3]/(fData.τ1[k]^2*sqrt(2)) + ν5[k]/(fData.τ1[k]^2*sqrt(2)) +
              cos(θδ[k])*vδ[k[2]]/fData.τ2[k]*(vmax[k[2]]/fData.τ2[k]*λlnc1[k] + vmin[k[2]]/fData.τ2[k]*λlnc2[k])/(fData.τ1[k]^2) for k in branchDict1[i])
              - sum(μ2[k,3]/(fData.τ2[k]^2*sqrt(2)) + ν2[k]/(fData.τ2[k]^2*sqrt(2)) + μ5[k,4]/(fData.τ2[k]^2*sqrt(2)) + ν5[k]/(fData.τ2[k]^2*sqrt(2)) +
              cos(θδ[k])*vδ[k[1]]/fData.τ1[k]*(vmax[k[1]]/fData.τ1[k]*λlnc1[k] + vmin[k[1]]/fData.τ1[k]*λlnc2[k])/(fData.τ2[k]^2) for k in branchDict2[i]) == 0);
  @constraint(sprob,vvConstr[k in fData.brList],
              λvv1[k]+λvv2[k]+λvv3[k]+λvv4[k] + λvve[k] - λvve[(k[2],k[1],k[3])] - μ2[k,1]
              - csmin[k]*λwc1[k] -csmax[k]*λwc2[k] - csmax[k]*λwc3[k] - csmin[k]*λwc4[k]
              - ssmin[k]*λws1[k]-ssmax[k]*λws2[k]-ssmax[k]*λws3[k]-ssmin[k]*λws4[k] == 0);
  @constraint(sprob,wcConstr[k in fData.brList],
              fData.g[k]*con1Par[k]*λptrans[k] - fData.b[k]*con1Par[k]*λqtrans[k] - μ5[k,1] +
              λwc1[k] + λwc2[k] + λwc3[k] + λwc4[k] + λwce[k] - λwce[(k[2],k[1],k[3])] -
              λtangent1[k]*tan(θDmax[k]) - λtangent2[k]*tan(θDmin[k]) +
              vδ[k[1]]*vδ[k[2]]/(fData.τ1[k]*fData.τ2[k])*cos(θϕ[k])*(λlnc1[k] + λlnc2[k]) == 0);
  @constraint(sprob,wsConstr[k in fData.brList],
              fData.b[k]*con1Par[k]*λptrans[k]+fData.g[k]*con1Par[k]*λqtrans[k] + λtangent1[k] + λtangent2[k] - μ5[k,2] +
              λws1[k] + λws2[k] + λws3[k] + λws4[k] + λwse[k] + λwse[(k[2],k[1],k[3])] +
              vδ[k[1]]*vδ[k[2]]/(fData.τ1[k]*fData.τ2[k])*sin(θϕ[k])*(λlnc1[k] + λlnc2[k]) == 0);
  @constraint(sprob,csConstr[k in fData.brList],
              λcs1[k] + λcs2[k] + λcs3[k] - μ3[k,2] + ν3[k] + λcse[k] - λcse[(k[2],k[1],k[3])]
              - vmin[k[1]]*vmin[k[2]]/(fData.τ1[k]*fData.τ2[k])*λwc1[k] - vmax[k[1]]*vmax[k[2]]*λwc2[k]/(fData.τ1[k]*fData.τ2[k])
              - vmin[k[1]]*vmin[k[2]]/(fData.τ1[k]*fData.τ2[k])*λwc3[k] - vmax[k[1]]*vmax[k[2]]*λwc4[k]/(fData.τ1[k]*fData.τ2[k]) == 0);
  @constraint(sprob,ssConstr[k in fData.brList],
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
  @constraint(sprob,θConstr[i in fData.IDList],
              λθRef*refBusInd[i] + sum(λθ1[k] + λθ2[k] - csConst[k]*λcs3[k] - cos(θu[k]/2)*λss4[k] - cos(θu[k]/2)*λss5[k]
              - sqrt(1-cos(θu[k]))/θu[k]*μ3[k,1] for k in branchDict1[i])
              + sum(-λθ1[k] - λθ2[k] + csConst[k]*λcs3[k] + cos(θu[k]/2)*λss4[k] + cos(θu[k]/2)*λss5[k]
              + sqrt(1-cos(θu[k]))/(θu[k])*μ3[k,1] for k in branchDict2[i]) == 0);
  @constraint(sprob,opConstr1[i in fData.IDList],λop1[i] >= λpi[i]);
  @constraint(sprob,opConstr2[i in fData.IDList],λop2[i] >= -λpi[i]);
  @constraint(sprob,oqConstr1[i in fData.IDList],λoq1[i] >= λqi[i]);
  @constraint(sprob,oqConstr2[i in fData.IDList],λoq2[i] >= -λqi[i]);

  # set up the objective function
  # rescale the objective function about really small GP0/GQ0, WPmax/WQmax
  conGPpar1 = Dict();
  conGQpar1 = Dict();
  conGPpar2 = Dict();
  conGQpar2 = Dict();
  for i in fData.IDList
    if (abs(uData[i].GP0) <= 1e-3)&(abs(uData[i].GP0) > 0)
      conGPpar1[i] = 1/sqrt(abs(uData[i].GP0));
      conGPpar2[i] = sqrt(abs(uData[i].GP0));
    else
      conGPpar1[i] = 1;
      conGPpar2[i] = 1;
    end
    if (abs(uData[i].GQ0) <= 1e-3)&(abs(uData[i].GQ0) > 0)
      conGQpar1[i] = 1/sqrt(abs(uData[i].GQ0));
      conGQpar2[i] = sqrt(abs(uData[i].GQ0));
    else
      conGQpar1[i] = 1;
      conGQpar2[i] = 1;
    end
  end
  @variable(sprob, λpis[i in fData.IDList]);
  @variable(sprob, λqis[i in fData.IDList]);
  @constraint(sprob,levelp[i in fData.IDList], λpis[i] == conGPpar2[i]*λpi[i]);
  @constraint(sprob,levelq[i in fData.IDList], λqis[i] == conGQpar2[i]*λqi[i]);

  M1 = 10;
  @variable(sprob, ropplus[i in fData.IDList] >= 0);
  @variable(sprob, ropminus[i in fData.IDList] >= 0);
  @variable(sprob, roqplus[i in fData.IDList] >= 0);
  @variable(sprob, roqminus[i in fData.IDList] >= 0);
  @constraint(sprob, ubRopplus1[i in fData.IDList], ropplus[i] <= ypplus[i]);
  @constraint(sprob, lbRopplus1[i in fData.IDList], ropplus[i] >= -ypplus[i]);
  @constraint(sprob, ubRopplus2[i in fData.IDList], ropplus[i] <= λop2[i] + (1-ypplus[i])*M1);
  @constraint(sprob, lbRopplus2[i in fData.IDList], ropplus[i] >= λop2[i] - (1-ypplus[i])*M1);
  @constraint(sprob, ubRopminus1[i in fData.IDList], ropminus[i] <= ypminus[i]);
  @constraint(sprob, lbRopminus1[i in fData.IDList], ropminus[i] >= -ypminus[i]);
  @constraint(sprob, ubRopminus2[i in fData.IDList], ropminus[i] <= λop2[i] + (1-ypminus[i])*M1);
  @constraint(sprob, lbRopminus2[i in fData.IDList], ropminus[i] >= λop2[i] - (1-ypminus[i])*M1);
  @constraint(sprob, ubRoqplus1[i in fData.IDList], roqplus[i] <= yqplus[i]);
  @constraint(sprob, lbRoqplus1[i in fData.IDList], roqplus[i] >= -yqplus[i]);
  @constraint(sprob, ubRoqplus2[i in fData.IDList], roqplus[i] <= λoq2[i] + (1-yqplus[i])*M1);
  @constraint(sprob, lbRoqplus2[i in fData.IDList], roqplus[i] >= λoq2[i] - (1-yqplus[i])*M1);
  @constraint(sprob, ubRoqminus1[i in fData.IDList], roqminus[i] <= yqminus[i]);
  @constraint(sprob, lbRoqminus1[i in fData.IDList], roqminus[i] >= -yqminus[i]);
  @constraint(sprob, ubRoqminus2[i in fData.IDList], roqminus[i] <= λoq2[i] + (1-yqminus[i])*M1);
  @constraint(sprob, lbRoqminus2[i in fData.IDList], roqminus[i] >= λoq2[i] - (1-yqminus[i])*M1);

  @expression(sprob,restTerm,-sum(fData.rateA[k]*ν1[k] for k in fData.brList if !(fData.rateA[k] == Inf)) -
              sum(rpplus[i]*(uData[i].GPmax - uData[i].GP0) + rpminus[i]*(-uData[i].GP0 + uData[i].GPmin) for i in fData.IDList if !(i in rmListP)) -
              sum(rqplus[i]*(uData[i].GQmax - uData[i].GQ0) + rqminus[i]*(-uData[i].GQ0 + uData[i].GQmin) for i in fData.IDList if !(i in rmListQ)) -
              sum(uData[i].GP0*conGPpar1[i]*λpis[i] + uData[i].GQ0*conGQpar1[i]*λqis[i] - μ4[i,2]/4 + ν4[i]/4 + λvu[i]*vmax[i] + λvl[i]*vmin[i] +
              λop1[i]*uData[i].WPmax + λoq1[i]*uData[i].WQmax - vmax[i]*vmin[i]*λv[i] + λop2[i]*(uData[i].WPmax + uData[i].RESP0) + λoq2[i]*(uData[i].WQmax + uData[i].RESQ0) +
              ropplus[i]*(uData[i].RESPmax - uData[i].RESP0) + ropminus[i]*(uData[i].RESPmin - uData[i].RESP0) +
              roqplus[i]*(uData[i].RESQmax - uData[i].RESQ0) + roqminus[i]*(uData[i].RESQmin - uData[i].RESQ0) for i in fData.IDList) -
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
  return sprob,restTerm;
end
