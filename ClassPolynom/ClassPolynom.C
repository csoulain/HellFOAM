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

#include "ClassPolynom.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ClassPolynom::ClassPolynom()
{}


Foam::ClassPolynom::ClassPolynom
(
//    const dictionary& polyProperties,
	const List<scalar> & polyCoeff,
    const word& polyName
)
:
	name_(polyName),
	polyCoeff_(polyCoeff)
//	polyCoeff_(polyProperties.lookup(polyName))
{}



Foam::ClassPolynom::ClassPolynom
(
    ClassPolynom & polyNome
)
{
	this->polyCoeff_ = polyNome.polyCoeff();
	this->name_  = polyNome.name();
}


Foam::ClassPolynom::ClassPolynom
(
    const ClassPolynom & polyNome
)
{
	this->polyCoeff_ = polyNome.polyCoeff();
	this->name_  = polyNome.name();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::ClassPolynom::~ClassPolynom()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::ClassPolynom::print()
{
	for(int k=0; k<polyCoeff_.size() ; k++)
	{
	Info << name_<<"["<<k<<"] = " << polyCoeff_[k] << endl;
	}		
	Info << nl << endl;
}

double Foam::ClassPolynom::poly( const double & X)
{
	double polyT = 0.0;			

	for(int i=0; i<polyCoeff_.size(); i++)
	{
		polyT += polyCoeff_[i]*Foam::pow(X,i);			
	}

	return polyT;
}

void Foam::ClassPolynom::poly( const volScalarField & X, volScalarField & polyX)
{
	polyX = 0.0*polyX;			

	forAll(polyX,cellI)
	{
		for(int i=0; i<polyCoeff_.size(); i++)
		{
			polyX[cellI] += polyCoeff_[i]*Foam::pow(X[cellI],i);			
		}
	}
	forAll(polyX.boundaryField(),cellI)
	{
		for(int i=0; i<polyCoeff_.size(); i++)
		{
			polyX.boundaryField()[cellI] += polyCoeff_[i]*Foam::pow(X.boundaryField()[cellI],i);			
		}
	}
}

Foam::ClassPolynom Foam::ClassPolynom::deriv()
{
	List<scalar> derivPolyCoeff(polyCoeff_.size()-1);

	derivPolyCoeff = 0;

	for(int k=0; k<derivPolyCoeff.size() ; k++)
	{
		derivPolyCoeff[k] = (k+1)*polyCoeff_[k+1];
	}

	ClassPolynom derivPoly( derivPolyCoeff,"derivPolyCoeff");  //"deriv"+name_;

	return derivPoly;
}


//Foam::ClassPolynom Foam::ClassPolynom::addPoly(const ClassPolynom& polyCoeff1, const ClassPolynom& polyCoeff2)
Foam::ClassPolynom Foam::ClassPolynom::addPoly( ClassPolynom& polyCoeff1, ClassPolynom& polyCoeff2)
{
	// attention, lorsque les degrés des Polynomes sont différents, l'algo suivant
	// est peut être source d'erreurs.
	// attention aussi lors que les termes de plus haut degré s'annulent.

	List<scalar> polyCoeffResult(max(polyCoeff1.polyCoeff().size(),polyCoeff2.polyCoeff().size()));
	
	for(int k=0; k<polyCoeffResult.size() ; k++)
	{
		polyCoeffResult[k]=polyCoeff1.polyCoeff()[k]+polyCoeff2.polyCoeff()[k];
	}

	ClassPolynom addPolyResult( polyCoeffResult, polyCoeff1.name()+"Plus"+polyCoeff2.name());

	return addPolyResult;
}

// ************************************************************************* //


