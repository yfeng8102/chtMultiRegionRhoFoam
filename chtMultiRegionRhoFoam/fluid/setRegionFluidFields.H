    //首先是变量名称的转化,将各个流体区域的变量转化为标准求解器中的形式
    fvMesh& mesh = fluidRegions[i];  //网格

    rhoThermo& thermo = thermoFluid[i];  //thermo
    thermo.validate(args.executable(), "h", "e");  //

    volScalarField& rho = rhoFluid[i];  //密度
    volVectorField& U = UFluid[i];  //速度
    surfaceScalarField& phi = phiFluid[i];  //通量

    compressible::turbulenceModel& turb = turbulence[i];  //湍流模型,注意与一般求解器中的turbulence不同,这里使用的是turb
    volScalarField& K = KFluid[i];  //动能
    volScalarField& dpdt = dpdtFluid[i];  //dpdt

    volScalarField& p = thermo.p();  //利用thermo计算压力
    const volScalarField& psi = thermo.psi();  //计算psi

    volScalarField& p_rgh = p_rghFluid[i];  //p_rgh

    const dimensionedVector& g = gFluid[i]; //重力
    const volScalarField& gh = ghFluid[i];  //gh
    const surfaceScalarField& ghf = ghfFluid[i];  //参考gh

    radiation::radiationModel& rad = radiation[i];  //辐射

    IOMRFZoneList& MRF = MRFfluid[i];  //MRF
    fv::options& fvOptions = fluidFvOptions[i];  //源项

    const dimensionedScalar initialMass  //初始质量
    (
        "initialMass",
        dimMass,
        initialMassFluid[i]
    );

    const bool frozenFlow = frozenFlowFluid[i];  //流体是否"冻结"

    bool& resReachedFluid = residualReachedFluid[i];  //控制参数
    bool& residualControlUsed = residualControlUsedFluid[i];  //控制参数
