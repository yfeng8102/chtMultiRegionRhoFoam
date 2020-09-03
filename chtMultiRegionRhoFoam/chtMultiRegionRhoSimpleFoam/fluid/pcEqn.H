{
    rho = thermo.rho();
    //the following is not in rhoSimpleFoam,
    //so deleted
	//restore to add rho limitor
    
    rho = max(rho, rhoMin[i]);
    rho = min(rho, rhoMax[i]);
	//according to the following, relax was commented
    //rho.relax();
    

    volScalarField rAU("rAU", 1.0/UEqn.A());
    //add rAtU
    volScalarField rAtU(1.0/(1.0/rAU - UEqn.H1()));
    //change to volScalarField
    //surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rho*rAU));
    volScalarField rhorAtU("rhorAtU", rho*rAtU);
    //change p_rgh to p
    //volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p_rgh));
    volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
    tUEqn.clear();

    //no buoyancy, so delete phig
    //surfaceScalarField phig(-rhorAUf*ghf*fvc::snGrad(rho)*mesh.magSf());

    //modify phiHbyA
    /*
    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        fvc::flux(rho*HbyA)
    );
    */
    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        fvc::interpolate(rho)*fvc::flux(HbyA)
    );

    MRF.makeRelative(fvc::interpolate(rho), phiHbyA);

    //change closedVlolume
    //bool closedVolume = adjustPhi(phiHbyA, U, p_rgh);
    bool closedVolume = false;

    //no buoyancy
    //phiHbyA += phig;

    // Update the pressure BCs to ensure flux consistency
    //change p_rgh to p
    //constrainPressure(p_rgh, rho, U, phiHbyA, rhorAUf, MRF);
    constrainPressure(p, rho, U, phiHbyA, rhorAtU, MRF);

    //rhoSimpleFoam considers compressibility, so delete
    /*
    dimensionedScalar compressibility = fvc::domainIntegrate(psi);
    bool compressible = (compressibility.value() > SMALL);
    */
    
    // Solve pressure
    //modify for transonic and subsonic fluid flow
    /*
    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix p_rghEqn
        (
            fvm::laplacian(rhorAUf, p_rgh) == fvc::div(phiHbyA)
        );

        p_rghEqn.setReference
        (
            pRefCell,
            compressible ? getRefCellValue(p_rgh, pRefCell) : pRefValue
        );

        solvPerfp_rgh = p_rghEqn.solve();

        if (nonOrth == nNonOrthCorr)
        {
            // Calculate the conservative fluxes
            phi = phiHbyA - p_rghEqn.flux();

            // Explicitly relax pressure for momentum corrector
            p_rgh.relax();

            // Correct the momentum source with the pressure gradient flux
            // calculated from the relaxed pressure
            U = HbyA + rAU*fvc::reconstruct((phig - p_rghEqn.flux())/rhorAUf);
            U.correctBoundaryConditions();
            fvOptions.correct(U);
        }
    }
    */
    if (transonic)
    {
        //add phid and modify HbyA as rhoSimpleFoam
        surfaceScalarField phid
        (
            "phid",
            (fvc::interpolate(psi)/fvc::interpolate(rho))*phiHbyA
        );

        phiHbyA +=
            fvc::interpolate(rho*(rAtU - rAU))*fvc::snGrad(p)*mesh.magSf()
          - fvc::interpolate(psi*p)*phiHbyA/fvc::interpolate(rho);
          
        HbyA -= (rAU - rAtU)*fvc::grad(p);
        
        //define and solve pEqn using original controls
        for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
        {
            fvScalarMatrix pEqn
            (
                fvc::div(phiHbyA)
              + fvm::div(phid, p)
              - fvm::laplacian(rhorAtU, p)
             ==
                fvOptions(psi, p, rho.name())
            );

            pEqn.relax();
            
            pEqn.setReference
            (
                pressureControl.refCell(),
                pressureControl.refValue()
            );
            
            solvPerfp = pEqn.solve();

            if (nonOrth == nNonOrthCorr)
            {
                phi = phiHbyA + pEqn.flux();
            }
        }
    }
    else
    {
        closedVolume = adjustPhi(phiHbyA, U, p);
        
        phiHbyA += fvc::interpolate(rho*(rAtU - rAU))*fvc::snGrad(p)*mesh.magSf();
        HbyA -= (rAU - rAtU)*fvc::grad(p);
        
        for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
        {
            fvScalarMatrix pEqn
            (
                fvc::div(phiHbyA)
              - fvm::laplacian(rhorAtU, p)
              ==
                fvOptions(psi, p, rho.name())
            );
            
            pEqn.setReference
            (
                pressureControl.refCell(),
                pressureControl.refValue()
            );
            
            solvPerfp = pEqn.solve();

            if (nonOrth == nNonOrthCorr)
            {
                phi = phiHbyA + pEqn.flux();
            }
        }
    }

    //no p_rgh
    //p = p_rgh + rho*gh;

    #include "continuityErrs.H"
    
    //move p relax and U correction here and add pLimited
    // Explicitly relax pressure for momentum corrector
    p.relax();

    U = HbyA - rAtU*fvc::grad(p);
    U.correctBoundaryConditions();
    fvOptions.correct(U);
    
    bool pLimited = pressureControl.limit(p);

    //modify 
    // For closed-volume cases adjust the pressure level
    // to obey overall mass continuity
    /*
    if (closedVolume && compressible)
    {
        p += (initialMass - fvc::domainIntegrate(thermo.rho()))
            /compressibility;
        p_rgh = p - rho*gh;
    }
    */
    if (closedVolume)
    {
        p += (initialMass - fvc::domainIntegrate(psi*p))
            /fvc::domainIntegrate(psi);
    }
    
    //add boundary correct
    if (pLimited || closedVolume)
    {
        p.correctBoundaryConditions();
    }
    

    rho = thermo.rho();
    //the following is deleted as rhoSimpleFoam
    //if this affect converging, then change it
	//restore
    rho = max(rho, rhoMin[i]);
    rho = min(rho, rhoMax[i]);

    //rho.relax();
	
	//change to if, transonic will not relax
	//accodding to this, relax should be commented
    if (!transonic)
    {
        rho.relax();
    }

    Info<< "Min/max rho:" << min(rho).value() << ' '
        << max(rho).value() << endl;
}
