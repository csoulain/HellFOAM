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

	Info << "----------------TEST DE LA CLASSE---------------- " << nl << endl;

        ClassPolynom polyRhon
        (
            polyProperties,
            "rhonCoeff"
        );
	
		Info << "Print rhon0 avec classe" << nl << endl;
		polyRhon.print();

		Info << "rhon avec classe = "<< polyRhon.poly(Tb.value()) << nl << endl;

	Info << "**************************************************" << nl << endl;

    List<scalar> rhonCoeff(polyProperties.lookup("rhonCoeff"));
	Info << "Print rhon" << nl << endl;
	printPoly(rhonCoeff);


    List<scalar> rhosCoeff(polyProperties.lookup("rhosCoeff"));
    List<scalar> rhoCoeff = addPoly(rhonCoeff,rhosCoeff);

    List<scalar>    sCoeff(polyProperties.lookup("sCoeff"));

    List<scalar>    kfCoeff(polyProperties.lookup("kCoeff"));

    List<scalar> muCoeff(polyProperties.lookup("muCoeff"));

    List<scalar> ACoeff(polyProperties.lookup("ACoeff"));

//	List<scalar> dSdTCoeff = derivPoly1(sCoeff);  

//   List<scalar> sRhoCoeff = multiplyPoly(rhoCoeff,sCoeff);

//   List<scalar> DsRhoDTCoeff = derivPoly1(sRhoCoeff);


    Info<< "Thermo update\n" << endl;

	s = calcPoly(sCoeff, Tb, s);
	Info << "s = " << s.value() << nl << endl;

	mun = calcPoly(muCoeff, Tb, mun);
	Info << "mun = " << mun.value() << nl << endl;

	rhon = calcPoly(rhonCoeff, Tb, rhon);
	Info << "rhon = " << rhon.value() << nl << endl;

	rhos = calcPoly(rhosCoeff, Tb, rhos);
	Info << "rhos = " << rhos.value() << nl << endl;

	Info << "rhos/rhon = " << rhos.value()/rhon.value() << nl << endl;

	rho.value() = rhos.value() + rhon.value();
	Info << "rho = " << rho.value() << nl << endl;

	kf = calcPoly(kfCoeff, Tb, kf);
	Info << "kf = " << kf.value() << nl << endl;

	A = calcPoly(ACoeff, Tb, A);
	Info << "A = " << A.value() << nl << endl;

//	dimensionedScalar dSdT ("dSdT", s.dimensions()/Tb.dimensions(), 0.0);
//	dSdT = calcPoly(dSdTCoeff, Tb, dSdT);
//	Info << "dSdT = " << dSdT.value() << nl << endl;

