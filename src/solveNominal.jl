# solve the nominal problem in convex and nonconvex settings
function solveNominalNC(fData,uData,vmax,vmin,θDmax,θDmin,βSwitch = 1)
    mp = ACOPF_NC(fData,uData,vmax,vmin,θDmax,θDmin,βSwitch);
    solve(mp);
    sphat = Dict();
    sqhat = Dict();
    for i in fData.genIDList
        sphat[i] = getvalue(mp[:sp][i]);
        sqhat[i] = getvalue(mp[:sq][i]);
    end
    cmp = getobjectivevalue(mp);

    return sphat,sqhat,cmp;
end

function solveNominalC(fData,uData,vmax,vmin,θDmax,θDmin,βSwitch = 1)
    mp = ACOPF_C_QP(fData,uData,vmax,vmin,θDmax,θDmin,βSwitch);
    solve(mp);
    sphat = Dict();
    sqhat = Dict();
    for i in fData.genIDList
        sphat[i] = getvalue(mp[:sp][i]);
        sqhat[i] = getvalue(mp[:sq][i]);
    end
    cmp = getobjectivevalue(mp);

    return sphat,sqhat,cmp;
end
