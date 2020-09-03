    // Solve the Momentum equation

    MRF.correctBoundaryVelocity(U);  //MRF

    tmp<fvVectorMatrix> tUEqn  //定义tUEqn方程
    (
        fvm::ddt(rho, U) + fvm::div(phi, U)
      + MRF.DDt(rho, U)
      + turb.divDevRhoReff(U)
     ==
        fvOptions(rho, U)
    );
    fvVectorMatrix& UEqn = tUEqn.ref();  //进而定义UEqn

    UEqn.relax();  //松弛
 
    fvOptions.constrain(UEqn);  //限制

    if (momentumPredictor)  //求解
    {
        solvPerfU =
            solve
            (
                UEqn
              ==
                fvc::reconstruct
                (
                    (
                      - ghf*fvc::snGrad(rho)
                      - fvc::snGrad(p_rgh)
                    )*mesh.magSf()
                ),
                mesh.solver(U.select(finalIter))
            );

        fvOptions.correct(U);  //修正源项
        K = 0.5*magSqr(U);  //计算动能
    }
