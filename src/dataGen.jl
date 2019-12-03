# create fData from the given file
function makeFData(fileAdd)
  # generate the uncertainty set
  busST,genST,brST,cST,uST,baseMVA = readMP(fileAdd);
  IDList,Vmax,Vmin,gs,bs,Vmag,Vang,Pd,Qd,bType = parsebusST(busST,baseMVA);
  genIDList,Loc,LocRev,Pmax,Pmin,Qmax,Qmin,Pg,Qg = parsegenST(genST,baseMVA);
  brList,brRev,g,b,bc,angmax,angmin,rateA,τ1,τ2,σ = parsebrST(brST,baseMVA);
  cp,cq = parsecST(cST,genIDList,baseMVA);
  fData = constFixed(baseMVA,bType,IDList,genIDList,brList,brRev,Vmax,Vmin,Loc,LocRev,Pmax,Pmin,Qmax,Qmin,
          gs,bs,Vmag,Vang,Pd,Qd,Pg,Qg,g,b,bc,angmax,angmin,rateA,τ1,τ2,σ,cp,cq);

  # specify the line power flow constraint if there is none
  θu = pi/3;
  for k in fData.brList
    if fData.rateA[k] == Inf
      fData.rateA[k] = sqrt(fData.g[k]^2+fData.b[k]^2)*max(fData.Vmax[k[1]],fData.Vmax[k[2]])*sqrt(fData.Vmax[k[1]]^2 + fData.Vmax[k[2]]^2 - 2*fData.Vmax[k[1]]*fData.Vmin[k[2]]*cos(θu));
    end
  end
  return fData;
end

function autouGen(fData,parameters)
  # Input: fData,baseMVA
  # Output: uData

  uData = Dict();

  # read in the bounds of recourse generation adjustment
  ωp = Dict();
  for i in fData.IDList
    if i in keys(parameters["ωp"])
      ωp[i] = parameters["ωp"][i];
    else
      ωp[i] = 0;
    end
  end
  ωq = Dict();
  for i in fData.IDList
    if i in keys(parameters["ωq"])
      ωq[i] = parameters["ωq"][i];
    else
      ωq[i] = 0;
    end
  end

  # generate the percentage perturbation of the demand
  αtop = Dict();
  αbot = Dict();
  if (typeof(parameters["αtop"]) == Float64)|(typeof(parameters["αtop"]) == Int64)
    for i in fData.IDList
      αtop[i] = parameters["αtop"];
    end
  else
    for i in fData.IDList
      if i in keys(parameters["αtop"])
        αtop[i] = parameters["αtop"][i];
      else
        αtop[i] = 0;
      end
    end
  end
  if (typeof(parameters["αbot"]) == Float64)|(typeof(parameters["αbot"]) == Int64)
    for i in fData.IDList
      αbot[i] = parameters["αbot"];
    end
  else
    for i in fData.IDList
      if i in keys(parameters["αbot"])
        αbot[i] = parameters["αbot"][i];
      else
        αbot[i] = 0;
      end
    end
  end

  DPmax = Dict();
  DPmin = Dict();
  DQmax = Dict();
  DQmin = Dict();
  for i in fData.IDList
    if fData.Pd[i] >= 0
      DPmax[i] = fData.Pd[i] * (1 + αtop[i]);
      DPmin[i] = fData.Pd[i] * (1 - αbot[i]);
    else
      DPmax[i] = fData.Pd[i] * (1 - αtop[i]);
      DPmin[i] = fData.Pd[i] * (1 + αbot[i]);
    end
    if fData.Qd[i] >= 0
      DQmax[i] = fData.Qd[i] * (1 + αtop[i]);
      DQmin[i] = fData.Qd[i] * (1 - αbot[i]);
    else
      DQmax[i] = fData.Qd[i] * (1 - αtop[i]);
      DQmin[i] = fData.Qd[i] * (1 + αbot[i]);
    end
  end

  # generate the uncertain supply information
  GPmax = Dict();
  GPmin = Dict();
  GP0 = Dict();
  GQmax = Dict();
  GQmin = Dict();
  GQ0 = Dict();
  for i in fData.IDList
    if i in keys(parameters["βunc"])
      GPmax[i] = (1 + parameters["βunc"][i][3])*parameters["βunc"][i][1];
      GPmin[i] = (1 - parameters["βunc"][i][4])*parameters["βunc"][i][1];
      GP0[i] = parameters["βunc"][i][1];
      GQmax[i] = (1 + parameters["βunc"][i][3])*parameters["βunc"][i][2];
      GQmin[i] = (1 - parameters["βunc"][i][4])*parameters["βunc"][i][2];
      GQ0[i] = parameters["βunc"][i][2];
    else
      GPmax[i] = 0;
      GPmin[i] = 0;
      GP0[i] = 0;
      GQmax[i] = 0;
      GQmin[i] = 0;
      GQ0[i] = 0;
    end
  end

    # append to uData
  for i in fData.IDList
    uData[i] = uncertainData((GPmax[i] - DPmin[i]),(GQmax[i] - DQmin[i]),(GPmin[i] - DPmax[i]),(GQmin[i] - DQmax[i]),
                  GP0[i]-fData.Pd[i],GQ0[i]-fData.Qd[i],ωp[i],ωq[i],GPmax[i],GQmax[i],GPmin[i],GQmin[i],GP0[i],GQ0[i],
                  1500,-1500,1500,-1500);
  end

  # return the generated uncertainty data
  return fData,uData;
end

function makeUAData(fData,ωp,ωq,αtop = 0,αbot = 0,genData = Dict())
  # generate the uncertainty set
  ωpDict = Dict();
  ωqDict = Dict();
  for i in fData.genIDList
    if fData.Loc[i] in keys(ωpDict)
      ωpDict[fData.Loc[i]] += ωp*(fData.Pmax[i] - fData.Pmin[i]);
      ωqDict[fData.Loc[i]] += ωq*(fData.Qmax[i] - fData.Qmin[i]);
    else
      ωpDict[fData.Loc[i]] = ωp*(fData.Pmax[i] - fData.Pmin[i]);
      ωqDict[fData.Loc[i]] = ωq*(fData.Qmax[i] - fData.Qmin[i]);
    end
  end
  paramsTest = Dict("αtop" => αtop, "αbot" => αbot, "ωp" => ωpDict, "ωq" => ωqDict, "βunc" => genData);
  fData,uData = autouGen(fData,paramsTest);
  return fData,uData;
end
