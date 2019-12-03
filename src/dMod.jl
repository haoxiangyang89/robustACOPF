# allow wiggle room in the second stage, including transformers

# this is the function to create the initial master problem
# remained unchanged from initMod.jl

function createM(fData::fixedData,uData::Dict{Any,Any},numericOpt = 0)
  # initialize the master problem
  # input: network parameters from the data variable
 mprob = Model(solver = GurobiSolver(NumericFocus = numericOpt, OutputFlag = 0, Threads = 20));

  # set up the variables: active/reactive power injection
  @variable(mprob, fData.Pmin[i] <= sp[i in fData.genIDList] <= fData.Pmax[i]);
  @variable(mprob, fData.Qmin[i] <= sq[i in fData.genIDList] <= fData.Qmax[i]);

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

  # set up the objective function
  @expression(mprob,objExpr,0);
  xp = Dict();
  xq = Dict();
  zp = Dict();
  zq = Dict();
  for i in fData.genIDList
    cpGen = fData.cp[i];
    if (cpGen.model == 1)
      if cpGen.n > 4
        # use anonymous constraints here since we don't really need to know the cost calculation procedure
        xp[i] = @variable(mprob,[1:(div(cpGen.n,2)-1)],Bin);
        zp[i] = @variable(mprob);
        @constraint(mprob,sp[i] <= cpGen.params[cpGen.n-1]);
        @constraint(mprob,sp[i] >= cpGen.params[1]);
        @constraint(mprob,sum(xp[i][j] for j in 1:(div(cpGen.n,2)-1)) == 1);
        @constraint(mprob,[j in 1:(div(cpGen.n,2)-1)],sp[i] <= (1-xp[i][j])*(cpGen.params[cpGen.n-1]-cpGen.params[1]+10)+cpGen.params[(j+1)*2]);
        @constraint(mprob,[j in 1:(div(cpGen.n,2)-1)],sp[i] >= (xp[i][j]-1)*(cpGen.params[cpGen.n-1]-cpGen.params[1]+10)+cpGen.params[j*2]);
        @constraint(mprob,[j in 1:(div(cpGen.n,2)-1)],zp[i] >= sp[i]*(cpGen.params[(j+1)*2] - cpGen.params[j*2])/
          (cpGen.params[(j+1)*2-1] - cpGen.params[j*2-1])+(cpGen.params[j*2]-(cpGen.params[(j+1)*2] - cpGen.params[j*2])/
          (cpGen.params[(j+1)*2-1] - cpGen.params[j*2-1])*cpGen.params[j*2-1]));
        objExpr += zp[j];
      elseif cpGen.n == 4
        @constraint(mprob,sp[i] <= cpGen.params[cpGen.n-1]);
        @constraint(mprob,sp[i] >= cpGen.params[1]);
        objExpr += (cpGen.params[4] - cpGen.params[2])/(cpGen.params[3] - cpGen.params[1])*sp[i] + cpGen.params[2];
      else
        error("Not enough parameters");
      end
    else
      for j in 1:cpGen.n
        objExpr += cpGen.params[j]*(sp[i]*fData.baseMVA)^(cpGen.n-j);
      end
    end
    if !(fData.cq == Dict())
      cqGen = fData.cq[i];
      if cqGen.model == 1
        if cqGen.n > 4
          xq[i] = @variable(mprob,[1:(div(cqGen.n,2)-1)],Bin);
          zq[i] = @variable(mprob);
          @constraint(mprob,sq[i] <= cqGen.params[cqGen.n-1]);
          @constraint(mprob,sq[i] >= cqGen.params[1]);
          @constraint(mprob,sum(xq[i][j] for j in 1:(div(cqGen.n,2)-1)) == 1);
          @constraint(mprob,[j in 1:(div(cqGen.n,2)-1)],sq[i] <= (1-xq[i][j])*(cqGen.params[cqGen.n-1]-cqGen.params[1]+10)+cqGen.qarams[(j+1)*2]);
          @constraint(mprob,[j in 1:(div(cqGen.n,2)-1)],sq[i] >= (xp[i][j]-1)*(cqGen.params[cqGen.n-1]-cqGen.params[1]+10)+cqGen.params[j*2]);
          @constraint(mprob,[j in 1:(div(cqGen.n,2)-1)],zq[i] >= sq[i]*(cqGen.params[(j+1)*2] - cqGen.params[j*2])/
            (cqGen.params[(j+1)*2-1] - cqGen.params[j*2-1])+(cqGen.params[j*2]-(cqGen.params[(j+1)*2] - cqGen.params[j*2])/
            (cqGen.params[(j+1)*2-1] - cqGen.params[j*2-1])*cqGen.params[j*2-1]));
          objExpr += zq[j];
        elseif cqGen.n == 4
          @constraint(mprob,sq[i] <= cqGen.params[cqGen.n-1]);
          @constraint(mprob,sq[i] >= cqGen.params[1]);
          objExpr += (cqGen.params[4] - cqGen.params[2])/(cqGen.params[3] - cqGen.params[1])*sq[i] + cqGen.params[2];
        else
          error("Not enough parameters");
        end
      else
        for j in 1:cqGen.n
          objExpr += cqGen.params[j]*(sq[i]*fData.baseMVA)^(cqGen.n-j);
        end
      end
    end
  end
  @objective(mprob,Min,objExpr);
  return mprob;
end

# we allow some wiggle room in the second stage, which will change the subproblem
function createS(fData::fixedData,uData::Dict{Any,Any},vmax,vmin,θDmax,θDmin,numericOpt = 0)
  # initialize the sub problem
  # obtian θu from the θDmax and θDmin
  θu = Dict();
  for k in fData.brList
    θu[k] = max(abs(θDmax[k]),abs(θDmin[k]));
  end
  # input: network parameters from the data variable
#  sprob = Model(solver = GurobiSolver(NumericFocus = 3,BarConvTol = 1e-6,MIPGap = 1e-6,BarQCPConvTol = 1e-6,
  # OptimalityTol = 1e-6,IntFeasTol = 1e-6,FeasibilityTol = 1e-6, OutputFlag = 0, Threads = 1));
  sprob = Model(solver = GurobiSolver(NumericFocus = numericOpt, OutputFlag = 0, Threads = 20));
  # sprob = Model(solver = MosekSolver(MSK_DPAR_MIO_TOL_REL_GAP = 1e-6,MSK_DPAR_INTPNT_TOL_PFEAS = 1e-6,
  #   MSK_DPAR_MIO_TOL_FEAS = 1e-6, MSK_IPAR_LOG = 0, MSK_IPAR_NUM_THREADS = 20));
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
  ## dual variables for cosine relaxation
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
  @variable(sprob,λvl[i in fData.IDList] >= 0);
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
  @variable(sprob,λoq1[i in fData.IDList] >= 0);
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
      # @constraint(sprob,yq[i] == 0);
      push!(rmListQ,i);
    end
  end
  @variable(sprob, yp[i in fData.IDList;!(i in rmListP)], Bin);
  @variable(sprob, rp[i in fData.IDList;!(i in rmListP)]);
  @variable(sprob, yq[i in fData.IDList;!(i in rmListQ)], Bin);
  @variable(sprob, rq[i in fData.IDList;!(i in rmListQ)]);

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
  @constraint(sprob, ubRp1[i in fData.IDList;!(i in rmListP)], rp[i] <= yp[i]);
  @constraint(sprob, lbRp1[i in fData.IDList;!(i in rmListP)], rp[i] >= -yp[i]);
  @constraint(sprob, ubRp2[i in fData.IDList;!(i in rmListP)], rp[i] <= λpi[i] + (1-yp[i])*M);
  @constraint(sprob, lbRp2[i in fData.IDList;!(i in rmListP)], rp[i] >= λpi[i] - (1-yp[i])*M);
  @constraint(sprob, ubRq1[i in fData.IDList;!(i in rmListQ)], rq[i] <= yq[i]);
  @constraint(sprob, lbRq1[i in fData.IDList;!(i in rmListQ)], rq[i] >= -yq[i]);
  @constraint(sprob, ubRq2[i in fData.IDList;!(i in rmListQ)], rq[i] <= λqi[i] + (1-yq[i])*M);
  @constraint(sprob, lbRq2[i in fData.IDList;!(i in rmListQ)], rq[i] >= λqi[i] - (1-yq[i])*M);

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
  @constraint(sprob,opConstr1[i in fData.IDList],λop1[i] >= -λpi[i]);
  @constraint(sprob,opConstr2[i in fData.IDList],λop1[i] >= λpi[i]);
  @constraint(sprob,oqConstr1[i in fData.IDList],λoq1[i] >= -λqi[i]);
  @constraint(sprob,oqConstr2[i in fData.IDList],λoq1[i] >= λqi[i]);

  # set up the objective function
  @expression(sprob,restTerm,-sum(fData.rateA[k]*ν1[k] for k in fData.brList if !(fData.rateA[k] == Inf)) -
            sum((uData[i].GPmax - uData[i].GPmin)*rp[i] for i in fData.IDList if !(i in rmListP)) -
            sum((uData[i].GQmax - uData[i].GQmin)*rq[i] for i in fData.IDList if !(i in rmListQ)) -
            sum(μ4[i,2]/4 + ν4[i]/4 + λvu[i]*vmax[i] - λvl[i]*vmin[i]
              +uData[i].GPmin*λpi[i] + uData[i].GQmin*λqi[i] + λop1[i]*uData[i].WPmax + λoq1[i]*uData[i].WQmax
              -vmax[i]*vmin[i]*λv[i] for i in fData.IDList) -
            sum(-3/4*μ3[k,2] + 5/4*ν3[k] + csmax[k]*λcs1[k] + csmin[k]*λcs2[k] + (cos(θDmin[k]) + csConst[k]*(fData.σ[k]*θDmin[k]))*λcs3[k]
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

function createSH(fData::fixedData,uData::Dict{Any,Any},yphat,yqhat,vmax,vmin,θDmax,θDmin)
  # initialize the sub problem
  # input: network parameters from the data variable
  # obtian θu from the θDmax and θDmin
  θu = Dict();
  for k in fData.brList
    θu[k] = max(abs(θDmax[k]),abs(θDmin[k]));
  end
#  sprob = Model(solver = GurobiSolver(NumericFocus = 3,BarConvTol = 1e-6,MIPGap = 1e-6,BarQCPConvTol = 1e-6,
  # OptimalityTol = 1e-6,IntFeasTol = 1e-6,FeasibilityTol = 1e-6, OutputFlag = 0, Threads = 20));
  sprob = Model(solver = GurobiSolver(NumericFocus = 3, OutputFlag = 0, Threads = 20));
  # sprob = Model(solver = MosekSolver(MSK_DPAR_MIO_TOL_REL_GAP = 1e-6,MSK_DPAR_INTPNT_TOL_PFEAS = 1e-6,
  #   MSK_DPAR_MIO_TOL_FEAS = 1e-6, MSK_IPAR_LOG = 0, MSK_IPAR_NUM_THREADS = 1));
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
  @variable(sprob,λoq1[i in fData.IDList] >= 0);
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
  @constraint(sprob,vhatConstr[i in fData.IDList], -μ4[i,2] - ν4[i] + λv[i]
              + λpi[i]*fData.gs[i] - λqi[i]*fData.bs[i]
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
  @constraint(sprob,opConstr1[i in fData.IDList],λop1[i] >= -λpi[i]);
  @constraint(sprob,opConstr2[i in fData.IDList],λop1[i] >= λpi[i]);
  @constraint(sprob,oqConstr1[i in fData.IDList],λoq1[i] >= -λqi[i]);
  @constraint(sprob,oqConstr2[i in fData.IDList],λoq1[i] >= λqi[i]);

  # set up the objective function
  @expression(sprob,restTerm,-sum(fData.rateA[k]*ν1[k] for k in fData.brList if !(fData.rateA[k] == Inf)) -
            sum((uData[i].GPmax - uData[i].GPmin)*λpi[i]*yphat[i]+(uData[i].GQmax - uData[i].GQmin)*λqi[i]*yqhat[i]
            +μ4[i,2]/4 + ν4[i]/4 + λvu[i]*vmax[i] - λvl[i]*vmin[i]
            +uData[i].GPmin*λpi[i] + uData[i].GQmin*λqi[i] + λop1[i]*uData[i].WPmax + λoq1[i]*uData[i].WQmax
            -vmax[i]*vmin[i]*λv[i] for i in fData.IDList) -
            sum(-3/4*μ3[k,2] + 5/4*ν3[k] + csmax[k]*λcs1[k] + csmin[k]*λcs2[k] + (cos(θDmin[k]) + csConst[k]*(fData.σ[k]*θDmin[k]))*λcs3[k]
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
