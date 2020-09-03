{
    volScalarField& he = thermo.he();

    volScalarField muEff("muEff", turb.muEff()); 
    volTensorField gradU("gradU", fvc::grad(U));
    volTensorField tau("tau", muEff*(dev2(gradU.T()) + gradU));

    fvScalarMatrix EEqn
    (
        fvm::div(phi, he)
      + (
            he.name() == "e"
          ? fvc::div(phi, volScalarField("Ekp", 0.5*magSqr(U) + p/rho))
          : fvc::div(phi, volScalarField("K", 0.5*magSqr(U)))
        )
      - fvm::laplacian(turb.alphaEff(), he)
      - fvc::div(tau & U)
     ==
    //no buoyancy and no radiation 
    //    rho*(U&g)
    //  + rad.Sh(thermo, he)
    //  + fvOptions(rho, he)
        fvOptions(rho, he)
    );

    EEqn.relax();

    fvOptions.constrain(EEqn);

    solvPerfE = EEqn.solve();

    fvOptions.correct(he);

    thermo.correct();
    //no radiation
    //rad.correct();

    Info<< "Min/max T:" << min(thermo.T()).value() << ' '
        << max(thermo.T()).value() << endl;
}
