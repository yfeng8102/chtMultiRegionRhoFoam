if (finalIter)
{
    mesh.data::add("finalIteration", true);
}

if (frozenFlow)  //如果"冻结",即流体静止,只求解能量方程
{
    #include "EEqn.H"
}
else  //否则,正常求解流动方程
{
    if (oCorr == 0)
    {
        #include "rhoEqn.H"  //密度方程
    }

    #include "UEqn.H"  //动量方程
    #include "EEqn.H"  //能量方程

    // --- PISO loop
    for (int corr=0; corr<nCorr; corr++)
    {
        #include "pEqn.H"  //压力修正方程
    }

    turb.correct();  //修正湍流模型边界

    rho = thermo.rho();  //更新密度
}

if (finalIter)
{
    mesh.data::remove("finalIteration");
}
