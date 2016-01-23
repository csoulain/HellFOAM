/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "HepakThermo.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


Foam::HepakThermo::HepakThermo
(
	const fvMesh& mesh,
	const dictionary& polyProperties,
	const volScalarField & T
)
:
	T_(T),
//	polyRhon_(polyProperties, "rhonCoeff"),
	polyRhon_(polyProperties.lookup("rhonCoeff"), "rhonCoeff"),
	rhon_
    (
        IOobject
        (
            "rhon",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
		dimensionedScalar("rhon",dimDensity,0.0)
    ),
//	polyRhos_(polyProperties, "rhosCoeff"),
	polyRhos_(polyProperties.lookup("rhosCoeff"), "rhosCoeff"),
	rhos_
    (
        IOobject
        (
            "rhos",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
		dimensionedScalar("rhos",dimDensity,0.0)
    ),
//	polyMun_(polyProperties, "muCoeff"),
	polyMun_(polyProperties.lookup("muCoeff"), "muCoeff"),
	Mun_
    (
        IOobject
        (
            "Mun",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
		dimensionedScalar("Mun",dimensionSet(1,-1,-1,0,0),0.0)
    ),
//	polyA_(polyProperties, "ACoeff"),
	polyA_(polyProperties.lookup("ACoeff"), "ACoeff"),
	A_
    (
        IOobject
        (
            "A",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
		dimensionedScalar("A",dimensionSet(-1,1,1,0,0),0.0)
    ),
//	polyS_(polyProperties, "sCoeff"),
	polyS_(polyProperties.lookup("sCoeff"), "sCoeff"),
	s_
    (
        IOobject
        (
            "s",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
		dimensionedScalar("s",dimensionSet(0,2,-2,-1,0),0.0)
    ),
//	polyK_(polyProperties, "kCoeff"),
	polyK_(polyProperties.lookup("kCoeff"), "kCoeff"),
	k_
    (
        IOobject
        (
            "k",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
		dimensionedScalar("k",dimensionSet(1,1,-3,-1,0),0.0)
    ),
	dSdT_
    (
        IOobject
        (
            "dSdT",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
//		dimensionedScalar("dSdT",dimensionSet(1,1,-3,-1,0),0.0)
		dimensionedScalar("dSdT",s_.dimensions()/T_.dimensions(),0.0)
    ),
	rho_
    (
        IOobject
        (
            "rho",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
		dimensionedScalar("rho",dimDensity,0.0)
    ),
	dRhodT_
    (
        IOobject
        (
            "dRhodT",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
		dimensionedScalar("dRhodT",rho_.dimensions()/T_.dimensions(),0.0)
    )
{
	polyRho_ = Foam::ClassPolynom::addPoly(polyRhon_,polyRhos_)	;
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


Foam::HepakThermo::~HepakThermo()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::HepakThermo::update()
{
	// calcule P(T)=a0 + a1*T + a2*T^2 + ....... + an*T^n 
	// où T est un champ scalair et (ai)i est le polynome

	polyRhon_.poly(T_,rhon_);
	polyRhos_.poly(T_,rhos_);
	polyMun_.poly(T_,Mun_);
	polyA_.poly(T_,A_);
	polyS_.poly(T_,s_);
	polyK_.poly(T_,k_);

	polyS_.deriv().poly(T_,dSdT_);

	polyRho_.poly(T_,rho_);

	polyRho_.deriv().poly(T_,dRhodT_);

}

// ************************************************************************* //
