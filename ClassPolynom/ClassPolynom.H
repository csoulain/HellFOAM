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

Class
    Foam::ClassPolynom

SourceFiles
    ClassPolynom.C

\*---------------------------------------------------------------------------*/

#ifndef ClassPolynom_H
#define ClassPolynom_H

#include "dictionary.H"
#include "dimensionedScalar.H"
//#include "fvPatchFields.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class phaseModel Declaration
\*---------------------------------------------------------------------------*/

class ClassPolynom
{

protected:
    // Private data

//        dictionary dict_;

        //- Name of phase
        word name_;

		List<scalar> polyCoeff_;

public:

    // Constructors
	ClassPolynom();

    ClassPolynom
    (
        //const dictionary& polyProperties,
   		const List<scalar> & polyCoeff,
        const word& polyName
    );

	ClassPolynom
	(
	    ClassPolynom &
	);

	ClassPolynom
	(
	    const ClassPolynom &
	);


    //- Destructor
    virtual ~ClassPolynom();


    // Member Functions

        const word& name() const
        {
			return name_;
        }

        const List<scalar>& polyCoeff() const
        {
			return polyCoeff_;
        }

		void print();


		// calcule P(T)=a0 + a1*T + a2*T^2 + ....... + an*T^n 
		// où T est un champ scalair et (ai)i est le polynome

		double poly( const double &);

//		Foam::fvPatchField<double> poly( const Foam::fvPatchField<double> &);

		void poly( const volScalarField & X, volScalarField & polyX);

//		ClassPolynom operator*(const ClassPolynom& polyCoeff1, const ClassPolynom& polyCoeff2);

		ClassPolynom deriv();

//		ClassPolynom addPoly(const ClassPolynom&, const ClassPolynom&);
		static ClassPolynom addPoly(ClassPolynom&, ClassPolynom&);
//		ClassPolynom& operator+(const ClassPolynom&, const ClassPolynom&);
};
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

