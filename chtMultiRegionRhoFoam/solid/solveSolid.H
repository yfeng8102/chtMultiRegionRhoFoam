//求解固体区域
if (finalIter)
{
    mesh.data::add("finalIteration", true);
}

{
    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix hEqn  //定义hEqn,这里的h是solidThermo中的h
        (
            fvm::ddt(betav*rho, h)
          - (
               thermo.isotropic()
             ? fvm::laplacian(betav*thermo.alpha(), h, "laplacian(alpha,h)")
             : fvm::laplacian(betav*taniAlpha(), h, "laplacian(alpha,h)")
            )
          ==
            fvOptions(rho, h)
        );

        hEqn.relax();  //松弛

        fvOptions.constrain(hEqn);  //限制

        solvPerfh = hEqn.solve(mesh.solver(h.select(finalIter)));  //求解

        fvOptions.correct(h);  //修正源项
    }
}

thermo.correct();

Info<< "Min/max T:" << min(thermo.T()).value() << ' '
    << max(thermo.T()).value() << endl;

if (finalIter)
{
    mesh.data::remove("finalIteration");
}
