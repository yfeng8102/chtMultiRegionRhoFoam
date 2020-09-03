// Initialise fluid field pointer lists
//首先是初始化流体区域场的指针
PtrList<rhoThermo> thermoFluid(fluidRegions.size());  //以thermo类为例,声明一个thermo类型的数组,数组长度为流体区域的个数,即每个流体区域的thermo
PtrList<volScalarField> rhoFluid(fluidRegions.size());
PtrList<volVectorField> UFluid(fluidRegions.size());
PtrList<surfaceScalarField> phiFluid(fluidRegions.size());
PtrList<uniformDimensionedVectorField> gFluid(fluidRegions.size());
PtrList<uniformDimensionedScalarField> hRefFluid(fluidRegions.size());
PtrList<volScalarField> ghFluid(fluidRegions.size());
PtrList<surfaceScalarField> ghfFluid(fluidRegions.size());
PtrList<compressible::turbulenceModel> turbulence(fluidRegions.size());
PtrList<volScalarField> p_rghFluid(fluidRegions.size());
PtrList<radiation::radiationModel> radiation(fluidRegions.size());
PtrList<volScalarField> KFluid(fluidRegions.size());
PtrList<volScalarField> dpdtFluid(fluidRegions.size());

List<scalar> initialMassFluid(fluidRegions.size());
List<bool> frozenFlowFluid(fluidRegions.size(), false);
List<bool> residualReachedFluid(fluidRegions.size(), true);
List<bool> residualControlUsedFluid(fluidRegions.size(), false);

PtrList<IOMRFZoneList> MRFfluid(fluidRegions.size());
PtrList<fv::options> fluidFvOptions(fluidRegions.size());

// Populate fluid field pointer lists
//下面为流体区域的指针赋值,也就是根据各个流体区域的设置文件设置相应的指针
forAll(fluidRegions, i)
{
    Info<< "*** Reading fluid mesh thermophysical properties for region "
        << fluidRegions[i].name() << nl << endl;
    //fluidRegions[i]表示网格,及通常求解器中的mesh
    Info<< "    Adding to thermoFluid\n" << endl;
    //设置thermo,这里使用的是rhoThermo
    thermoFluid.set
    (
        i,
        rhoThermo::New(fluidRegions[i]).ptr()
    );

    Info<< "    Adding to rhoFluid\n" << endl;
    rhoFluid.set  //设置流体密度
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "rho",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            thermoFluid[i].rho()  //密度是由thermo模型得到的
        )
    );

    Info<< "    Adding to UFluid\n" << endl;
    UFluid.set  //设置流体速度
    (
        i,
        new volVectorField
        (
            IOobject
            (
                "U",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i]  //由网格读出
        )
    );

    Info<< "    Adding to phiFluid\n" << endl;
    phiFluid.set  //设置流体的通量phi
    (
        i,
        new surfaceScalarField
        (
            IOobject
            (
                "phi",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            linearInterpolate(rhoFluid[i]*UFluid[i])
                & fluidRegions[i].Sf()  //是计算得到的,rho * U \cdot Sf
        )
    );

    Info<< "    Adding to gFluid\n" << endl;
    gFluid.set  //设置流体的重力场,应该是通过文件读取,和重力项有关
    (
        i,
        new uniformDimensionedVectorField
        (
            IOobject
            (
                "g",
                runTime.constant(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        )
    );

    Info<< "    Adding to hRefFluid\n" << endl;
    hRefFluid.set  //设置hRef,也就是参考高度,和重力项有关
    (
        i,
        new uniformDimensionedScalarField
        (
            IOobject
            (
                "hRef",
                runTime.constant(),
                fluidRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            dimensionedScalar("hRef", dimLength, 0)
        )
    );

    dimensionedScalar ghRef  //计算ghRef,也就是参考高度处的重力势能,和重力项有关
    (
        mag(gFluid[i].value()) > SMALL
      ? gFluid[i]
          & (cmptMag(gFluid[i].value())/mag(gFluid[i].value()))*hRefFluid[i]
      : dimensionedScalar("ghRef", gFluid[i].dimensions()*dimLength, 0)
    );  //如果重力大于某个定值,则有ghRef = gh,否则有ghRef = 0(?)

    Info<< "    Adding to ghFluid\n" << endl;
    ghFluid.set  //设置gh场,和重力项有关,这里是网格中心的值
    (
        i,
        new volScalarField  //由定义可知定义的gh是网格中心处的值
        (
            "gh",
            (gFluid[i] & fluidRegions[i].C()) - ghRef  //gh相对值 = 绝对值 - 参考值
        )  //这样看来,fluidRegion[i].C()表示网格中心处的值
    );

    Info<< "    Adding to ghfFluid\n" << endl;
    ghfFluid.set  //设置gh场,和重力项有关,这里是网格面上的值
    (
        i,
        new surfaceScalarField  //由定义可知定义的ghf是网格面上的值
        (
            "ghf",
            (gFluid[i] & fluidRegions[i].Cf()) - ghRef
        )  //这样看来,fluidRegions[i].Cf())表示的是网格面上的值
    );

    Info<< "    Adding to turbulence\n" << endl;
    turbulence.set  //定义湍流模型
    (
        i,
        compressible::turbulenceModel::New
        (
            rhoFluid[i],
            UFluid[i],
            phiFluid[i],
            thermoFluid[i]
        ).ptr()
    );

    p_rghFluid.set  //定义p_rgh
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "p_rgh",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i]  //由网格读取
        )
    );

    // Force p_rgh to be consistent with p
    p_rghFluid[i] = thermoFluid[i].p() - rhoFluid[i]*ghFluid[i];  //修正p_rgh,使之与p协调

    fluidRegions[i].setFluxRequired(p_rghFluid[i].name());

    radiation.set  //设置辐射项
    (
        i,
        radiation::radiationModel::New(thermoFluid[i].T())
    );

    initialMassFluid[i] = fvc::domainIntegrate(rhoFluid[i]).value();

    Info<< "    Adding to KFluid\n" << endl;
    KFluid.set  //设置动能
    (
        i,
        new volScalarField  //定义在网格中心
        (
            "K",
            0.5*magSqr(UFluid[i])  //动能K = 0.5 * U^2
        )
    );

    Info<< "    Adding to dpdtFluid\n" << endl;
    dpdtFluid.set  //定义dpdt,求解h方程时有这一项
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "dpdt",
                runTime.timeName(),
                fluidRegions[i]
            ),
            fluidRegions[i],
            dimensionedScalar
            (
                "dpdt",
                thermoFluid[i].p().dimensions()/dimTime,
                0
            )
        )
    );
    //下面是pinple算法的一些设置
    const dictionary& pimpleDict =
        fluidRegions[i].solutionDict().subDict("PIMPLE");
    pimpleDict.readIfPresent("frozenFlow", frozenFlowFluid[i]);

    Info<< "    Adding MRF\n" << endl;
    MRFfluid.set  //设置MRF
    (
        i,
        new IOMRFZoneList(fluidRegions[i])
    );

    Info<< "    Adding fvOptions\n" << endl;
    fluidFvOptions.set  //设置FvOpetions
    (
        i,
        new fv::options(fluidRegions[i])
    );

    turbulence[i].validate();
}
