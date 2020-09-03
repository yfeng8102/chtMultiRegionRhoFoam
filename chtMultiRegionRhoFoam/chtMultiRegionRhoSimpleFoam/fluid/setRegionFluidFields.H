    const fvMesh& mesh = fluidRegions[i];

    //rhoThermo& thermo = thermoFluid[i];
    fluidThermo& thermo = thermoFluid[i];
    thermo.validate(args.executable(), "h", "e");

    volScalarField& rho = rhoFluid[i];
    volVectorField& U = UFluid[i];
    surfaceScalarField& phi = phiFluid[i];

    compressible::turbulenceModel& turb = turbulence[i];

    //here p is defined again, but this will not
    //affect pFluid[i]
    volScalarField& p = thermo.p();
    const volScalarField& psi = thermo.psi();

    IOMRFZoneList& MRF = MRFfluid[i];
    fv::options& fvOptions = fluidFvOptions[i];

    //pressureControl
    //define simpleDict again
    const dictionary& simpleDict =
        fluidRegions[i].solutionDict().subDict("SIMPLE");
    pressureControl pressureControl(p, rho, simpleDict);
    
    const dimensionedScalar initialMass
    (
        "initialMass",
        dimMass,
        initialMassFluid[i]
    );

    //no radiation, so deleted
    //radiation::radiationModel& rad = radiation[i];

    //p ref is deleted
    //const label pRefCell = pRefCellFluid[i];
    //const scalar pRefValue = pRefValueFluid[i];
    const bool frozenFlow = frozenFlowFluid[i];

    //p_rgh is deleted
    //volScalarField& p_rgh = p_rghFluid[i];

    //buoyancy variables are deleted
    //const dimensionedVector& g = gFluid[i];
    //const volScalarField& gh = ghFluid[i];
    //const surfaceScalarField& ghf = ghfFluid[i];

    bool& resReachedFluid = residualReachedFluid[i];
    bool& residualControlUsed = residualControlUsedFluid[i];
