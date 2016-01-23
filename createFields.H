    Info<< "Reading transportProperties\n" << endl;

    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

	bool activeGorterMellink
	(
		readBool(transportProperties.lookup("activeGorterMellink"))
	);
	scalar aGM = 1.0;
	if(activeGorterMellink == false)
	{
		aGM = 0.0;
	}

	bool activeGradW
	(
		readBool(transportProperties.lookup("activeGradW"))
	);

	scalar aGW = 1.0;
	if (activeGradW == false)
	{
		aGW = 0.0;
	}

    dimensionedScalar tau
    (
        transportProperties.lookup("tau")
    );

    dimensionedScalar mus
    (
        transportProperties.lookup("mus")
    );

    dimensionedScalar Q0
    (
        transportProperties.lookup("Q0")
    );

    dimensionedScalar Q ("Q", 0.0*Q0);


    Info<< "Reading field alpha\n" << endl;
    volScalarField alpha
    (
        IOobject
        (
            "alpha",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
	dimensionedScalar Valpha ("Valpha", fvc::domainIntegrate(alpha));

    Info<< "Reading field gamma\n" << endl;
    volScalarField gamma
    (
        IOobject
        (
            "gamma",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
		dimensionedScalar("gamma",dimensionSet(1,-3,-1,0,0,0,0),0.0)
    );

    Info<< "Reading field continuity\n" << endl;
    volScalarField continuity
    (
        IOobject
        (
            "continuity",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
		dimensionedScalar("continuity",dimensionSet(1,-3,-1,0,0,0,0),0.0)
    );


    Info<< "Reading field T\n" << endl;
    volScalarField T
    (
        IOobject
        (
            "T",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

/*
	volScalarField T0 ("T0",T);

    Info<< "Reading field dT\n" << endl;
    volScalarField dT
    (
        IOobject
        (
            "dT",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (T-T0)*1e9
    );
*/
    IOdictionary polyProperties
    (
        IOobject
        (
            "polyProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

	HepakThermo thermo(mesh, polyProperties, T);

	volScalarField& rhon = thermo.rhon();

	volScalarField& rhos = thermo.rhos();

	volScalarField& s = thermo.s();

	volScalarField& kf = thermo.kf();

	volScalarField& A = thermo.A();

	volScalarField& mun = thermo.mun();


	volScalarField& dSdT = thermo.dSdT();

	volScalarField& rho = thermo.rho();

	volScalarField& dRhodT = thermo.dRhodT();

	thermo.update();

    Info<< "Reading field p\n" << endl;
    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field Un\n" << endl;
    volVectorField Un
    (
        IOobject
        (
            "Un",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

	surfaceScalarField phin ("phin", fvc::interpolate(rhon)*linearInterpolate(Un) & mesh.Sf());

    Info<< "Reading field Us\n" << endl;
    volVectorField Us
    (
        IOobject
        (
            "Us",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

	surfaceScalarField phis ("phis", fvc::interpolate(rhos)*linearInterpolate(Us) & mesh.Sf());

    Info<< "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        (rhos*Us + rhon*Un)/rho
    );

    Info<< "Reading field phi\n" << endl;
    surfaceScalarField phi
    (
        IOobject
        (
            "phi",
            runTime.timeName(),
            mesh
        ),
//        fvc::interpolate(rhon/rho)*phin + fvc::interpolate(rhos/rho)*phis
		phin+phis
    );


    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p, mesh.solutionDict().subDict("PIMPLE"), pRefCell, pRefValue);

