    const dictionary& simple = fluidRegions[i].solutionDict().subDict("SIMPLE");

    const int nNonOrthCorr =
        simple.lookupOrDefault<int>("nNonOrthogonalCorrectors", 0);
    
    //there is no if for momentumPredictor,
    //so it is always true
    const bool momentumPredictor =
        simple.lookupOrDefault("momentumPredictor", true);
    
    
    //control transonic
    const bool transonic =
        simple.lookupOrDefault("transonic", false);
        
    //control consistent
    const bool consistent =
        simple.lookupOrDefault("consistent", false);
        
    //print SIMPLEControls
    /*
    Info<< "*** print SIMPLEControls: " << nl << endl;
    Info<< "nNonOrthogonalCorrectors: " << nNonOrthCorr << nl << endl;
    Info<< "momentumPredictor: " << momentumPredictor << nl << endl;
    Info<< "transonic: " << transonic << nl << endl;
    Info<< "consistent: " << consistent << nl << endl;
    */
  
        
    