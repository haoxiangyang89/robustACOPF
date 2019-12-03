# obtained the largest allowed α with the given ωp and ωq under 1-norm uncertainty set
function limitAlphaNC(fData,ωp,ωq,testMode,vmax,vmin,θDmax,θDmin,partCorr,KK = 1,ϵ = 1e-4,βunc = Dict())
  # build the dictionary of ωp and ωq
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
  # read in the uncertainty data and parse them
  feasBool = true;
  dList = [];
  ii = 0;
  parameters = Dict("αtop" => ii, "αbot" => KK*ii, "ωp" => ωpDict, "ωq" => ωqDict, "βunc" => βunc);
  fData,uData = autouGen(fData,parameters);
  Ω,ypo,yqo = genPartialG(fData,uData,partCorr,testMode);

  # calculate the nonconvex robust cost from the extensive formulation
  mext = solveExtensiveNConv(fData,uData,Ω,ypo,yqo,vmax,vmin,θDmax,θDmin);
  meStatus = solve(mext);
  if meStatus != :Optimal
    ncC = Inf;
    sphat = Dict();
    sqhat = Dict();
  else
    ncC = getobjectivevalue(mext);
    sphat = Dict();
    sqhat = Dict();
    for i in fData.genIDList
      sphat[i] = getvalue(mext[:sp][i]);
      sqhat[i] = getvalue(mext[:sq][i]);
    end
  end
  push!(dList,[ii,sphat,sqhat,ncC]);

  iil = 0;
  iiu = 1;
  ii = (iil + iiu)/2;
  while feasBool
    parameters = Dict("αtop" => ii, "αbot" => KK*ii, "ωp" => ωpDict, "ωq" => ωqDict, "βunc" => βunc);
    fData,uData = autouGen(fData,parameters);

    # calculate the nonconvex robust cost from the extensive formulation
    mext = solveExtensiveNConv(fData,uData,Ω,ypo,yqo,vmax,vmin,θDmax,θDmin);
    meStatus = solve(mext);
    if meStatus != :Optimal
      ncC = Inf;
      sphat = Dict();
      sqhat = Dict();
      iiu = ii;
      ii = (iil + iiu)/2;
    else
      ncC = getobjectivevalue(mext);
      sphat = Dict();
      sqhat = Dict();
      for i in fData.genIDList
        sphat[i] = getvalue(mext[:sp][i]);
        sqhat[i] = getvalue(mext[:sq][i]);
      end
      push!(dList,[ii,sphat,sqhat,ncC]);
      println(ii);
      iil = ii;
      ii = (iil + iiu)/2;
    end

    if iiu - iil < ϵ
      feasBool = false;
    end
  end
  return dList,iil;
end

# obtained the largest allowed α with the given ωp and ωq under testMode-norm uncertainty set
function limitAlphaC(fData,ωp,ωq,testMode,vmax,vmin,θDmax,θDmin,partCorr,KK = 1,ϵ = 1e-4,βunc = Dict())
  # build the dictionary of ωp and ωq
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
  # read in the uncertainty data and parse them
  feasBool = true;
  dList = [];
  ii = 0;
  parameters = Dict("αtop" => ii, "αbot" => KK*ii, "ωp" => ωpDict, "ωq" => ωqDict, "βunc" => βunc);
  fData,uData = autouGen(fData,parameters);

  try
    sphat,sqhat,nC,C,vio,lbCost = solveProcess_Budget_CorrB(fData,uData,vmax,vmin,θDmax,θDmin,partCorr,testMode,1e-4,2,0);
    push!(dList,[ii,sphat,sqhat,lbCost]);
  catch
    lbCost = Inf;
    sphat = Dict();
    sqhat = Dict();
    push!(dList,[ii,sphat,sqhat,lbCost]);
  end

  iil = 0;
  iiu = 1;
  ii = (iil + iiu)/2;
  while feasBool
    parameters = Dict("αtop" => ii, "αbot" => KK*ii, "ωp" => ωpDict, "ωq" => ωqDict, "βunc" => βunc);
    fData,uData = autouGen(fData,parameters);
    try
      sphat,sqhat,nC,C,vio,lbCost = solveProcess_Budget_CorrB(fData,uData,vmax,vmin,θDmax,θDmin,partCorr,testMode,1e-4,2,0);
      push!(dList,[ii,sphat,sqhat,lbCost]);
      println(ii);
      iil = ii;
      ii = (iil + iiu)/2;
    catch
      iiu = ii;
      ii = (iil + iiu)/2;
    end
    if iiu - iil < ϵ
      feasBool = false;
    end
  end
  return dList,iil;
end

# obtain the αmax for the convex and nonconvex setting
function getAlphaNC(fData,βList,testMode,vmax,vmin,θDmax,θDmin,partCorr,KK = 1,ϵ = 1e-4,βunc = Dict())
    αmaxDict = Dict();
    recordDict = Dict();
    for β in βList
        ωp = β;
        ωq = β;
        recordDict[β],αmaxDict[β] = limitAlphaNC(fData,ωp,ωq,testMode,vmax,vmin,θDmax,θDmin,partCorr,KK,ϵ,βunc);
        println("$(β) nonconvex is completed.");
    end
    return recordDict,αmaxDict;
end

function getAlphaC(fData,βList,testMode,vmax,vmin,θDmax,θDmin,partCorr,KK = 1,ϵ = 1e-4,βunc = Dict())
    αmaxDict = Dict();
    recordDict = Dict();
    for β in βList
        ωp = β;
        ωq = β;
        recordDict[β],αmaxDict[β] = limitAlphaC(fData,ωp,ωq,testMode,vmax,vmin,θDmax,θDmin,partCorr,KK,ϵ,βunc);
        println("$(β) convex is completed.");
    end
    return recordDict,αmaxDict;
end
