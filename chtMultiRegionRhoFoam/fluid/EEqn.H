//能量方程
{
    volScalarField& he = thermo.he();
    
    fvScalarMatrix EEqn  //定义EEqn,可以使用h或者e
    (
        fvm::ddt(rho, he) + fvm::div(phi, he)
      + fvc::ddt(rho, K) + fvc::div(phi, K)
      + (
            he.name() == "e"
          ? fvc::div
            (
                fvc::absolute(phi/fvc::interpolate(rho), U),
                p,
                "div(phiv,p)"
            )
          : -dpdt
        )
      - fvm::laplacian(turb.alphaEff(), he)
     ==
        rho*(U&g)
      + rad.Sh(thermo, he)
      + fvOptions(rho, he)
    );

    EEqn.relax();  //松弛

    fvOptions.constrain(EEqn);  //限制

    solvPerfE = EEqn.solve(mesh.solver(he.select(finalIter)));  //求解

    fvOptions.correct(he);  //修正源项

    thermo.correct();
    rad.correct();

    Info<< "Min/max T:" << min(thermo.T()).value() << ' '
        << max(thermo.T()).value() << endl;
}
