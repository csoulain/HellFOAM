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
    Foam::HepakThermo

SourceFiles
    HepakThermo.C

\*---------------------------------------------------------------------------*/

#ifndef HepakThermo_H
#define HepakThermo_H

#include "dictionary.H"
#include "volFields.H"
#include "dimensionedScalar.H"
#include "../ClassPolynom/ClassPolynom.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class HepakThermo Declaration
\*---------------------------------------------------------------------------*/

class HepakThermo
{
	protected:

//		List<scalar> rhonCoeff_;
		const volScalarField & T_;

		ClassPolynom polyRhon_;
		
		volScalarField rhon_;

		ClassPolynom polyRhos_;
		
		volScalarField rhos_;

		ClassPolynom polyMun_;
		
		volScalarField Mun_;

		ClassPolynom polyA_;
		
		volScalarField A_;

		ClassPolynom polyS_;
		
		volScalarField s_;

		ClassPolynom polyK_;
		
		volScalarField k_;

		volScalarField dSdT_;

		ClassPolynom polyRho_;

		volScalarField rho_;

		volScalarField dRhodT_;

	public:
	
		// constructor
		HepakThermo
		(
            		const fvMesh& mesh,
			const dictionary& polyProperties,
			const volScalarField & T
		);

		// destructor
		~HepakThermo();

		// function member

        volScalarField& rhon()
        {
            return rhon_;
        }

        volScalarField& rhos()
        {
            return rhos_;
        }

        volScalarField& mun()
        {
            return Mun_;
        }

        volScalarField& A()
        {
            return A_;
        }

        volScalarField& s()
        {
            return s_;
        }

        volScalarField& kf()
        {
            return k_;
        }


        volScalarField& dSdT()
        {
            return dSdT_;
        }


        volScalarField& rho()
        {
            return rho_;
        }

        volScalarField& dRhodT()
        {
            return dRhodT_;
        }


		void update();

};
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

