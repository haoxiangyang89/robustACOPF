# tighten the formulation in the first stage scenario free

function processTighten_SFP(fData,uData,ϵ = 1e-2,iMax = 20)
  # for each scenario
  vmax = Dict();
  vmin = Dict();
  θDmax = Dict();
  θDmin = Dict();
  for i in fData.IDList
    vmax[i] = fData.Vmax[i];
    vmin[i] = fData.Vmin[i];
  end
  for k in fData.brList
    θDmax[k] = pi/6;
    θDmin[k] = -pi/6;
  end

  stopBool = false;
  counter = 0;
  # tighten the bound until it converges to one value
  while !stopBool
    counter += 1;
    println("--- Iteration $(counter) ---");
    # record the bound from last iteration
    @everywhere vmax0 = Dict();
    @everywhere vmin0 = Dict();
    @everywhere θDmax0 = Dict();
    @everywhere θDmin0 = Dict();
    for i in fData.IDList
      vmax0[i] = vmax[i];
      vmin0[i] = vmin[i];
    end
    for k in fData.brList
      θDmax0[k] = θDmax[k];
      θDmin0[k] = θDmin[k];
    end
    # for each bus
    idList = SharedArray{Int64}((length(fData.IDList)));
    for i in 1:length(fData.IDList)
       idList[i] = fData.IDList[i];
    end
    vmaxTemp = pmap(i -> solveBV_SF(fData,uData,idList[i],1,vmax0,vmin0,θDmax0,θDmin0), 1:length(fData.IDList));
    vminTemp = pmap(i -> solveBV_SF(fData,uData,idList[i],0,vmax0,vmin0,θDmax0,θDmin0), 1:length(fData.IDList));
    for i in 1:length(fData.IDList)
      vmax[fData.IDList[i]] = vmaxTemp[i];
      vmin[fData.IDList[i]] = vminTemp[i];
    end
    # for each line
    brList = SharedArray{Int64,2}((length(fData.brList),3));
    for k in 1:length(fData.brList)
      brList[k,1] = fData.brList[k][1];
      brList[k,2] = fData.brList[k][2];
      brList[k,3] = fData.brList[k][3];
    end
    θDmaxTemp = pmap(k -> solveBtheta_SF(fData,uData,(brList[k,1],brList[k,2],brList[k,3]),1,vmax0,vmin0,θDmax0,θDmin0), 1:length(fData.brList));
    θDminTemp = pmap(k -> solveBtheta_SF(fData,uData,(brList[k,1],brList[k,2],brList[k,3]),0,vmax0,vmin0,θDmax0,θDmin0), 1:length(fData.brList));
    for k in 1:length(fData.brList)
      θDmax[fData.brList[k]] = θDmaxTemp[k];
      θDmin[fData.brList[k]] = θDminTemp[k];
    end
    if counter < iMax
      stopBool = true;
      for i in fData.IDList
        if abs(vmax0[i] - vmax[i]) > ϵ
          stopBool = false;
        end
        if abs(vmin0[i] - vmin[i]) > ϵ
          stopBool = false;
        end
      end
      for k in fData.brList
        if abs(θDmax0[k] - θDmax[k]) > ϵ
          stopBool = false;
        end
        if abs(θDmin0[k] - θDmin[k]) > ϵ
          stopBool = false;
        end
      end
    else
      stopBool = true;
    end
  end

  for i in fData.IDList
    if vmax[i] - vmin[i] <= 1e-6
      vmax[i] += 1e-5;
      vmin[i] -= 1e-5;
    end
  end

  for k in fData.brList
    if θDmax[k] - θDmin[k] <= 1e-6
      θDmax[k] += 1e-5;
      θDmin[k] -= 1e-5;
    end
  end

  for k in fData.brList
    if (abs(θDmax[k] + θDmin[k]) > 0)&(abs(θDmax[k] + θDmin[k]) < 1e-3)
      if abs(θDmax[k]) < abs(θDmin[k])
        θDmax[k] = abs(θDmin[k]);
      else
        θDmin[k] = -abs(θDmax[k]);
      end
    end
  end
  return vmax,vmin,θDmax,θDmin;
end
