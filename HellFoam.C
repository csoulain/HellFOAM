/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2011 OpenCFD Ltd.
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

Application
    HellFoam

Description
    Solve the He(II) superfluid equations. The solver is based on the Super-PISO
	algorithm. 

    If you use this program, please cite
	"
	Soulaine, C.; Quintard, M.; Allain, H.; Baudouy, B. & Van Weelderen, R. 
	A PISO-like algorithm to simulate superfluid helium flow with the two-fluid model 
	Computer Physics Communications , 2015, 187, 20-28
	"

Authors
	2013-07-25: Cyprien Soulaine (CS, cyprien.soulaine@gmail.com) - First version
	2013-09-09: CS - Second major release
					* implementation of temperature-dependant variable. 
					* Correction of the temperature equation.
					* Upgrade to OpenFOAM 2.2.1
	2014-03-24: CS - minor release
					* Upgrade to OpenFOAM 2.3.0 
	2016-01-08: CS - minor release
					* Upgrade to OpenFOAM 3.0.1 

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

#include "ClassPolynom/ClassPolynom.H"
#include "HepakThermo/HepakThermo.H"
#include "pimpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"

    #include "createTime.H"
    #include "createMesh.H"
//    #include "readGravitationalAcceleration.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"
    #include "CourantNo.H"

    pimpleControl pimple(mesh);


    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {

        #include "CourantNos.H"

        runTime++;
        Info<< "Time = " << runTime.timeName() << nl << endl;

		thermo.update();

		Q = Q0*(scalar(1)-Foam::exp(-runTime.value()/tau.value()));

//		Q = Q0;

#if 0
		if ((runTime.value() <= scalar(20.e-3)) && (runTime.value() >= scalar(0.1e-3)) ) 
		{
	//		Q = Q0*(scalar(1)-Foam::exp(-runTime.value()/tau.value())); // à tester : sans rampe
			Q = Q0;
		}		
		else
		{
			Q = 0.0*Q;
		}

		//Q = Q0;
#endif

		Info << " max(T) = " << max(T).value() << " K" << nl << endl;
		Info << " Puissance de flux (watt) = " << Q.value() <<" sur " << Q0.value() << " W" << nl << endl;


        // --- Pressure-velocity PIMPLE corrector loop
		while (pimple.loop())
        {

			volScalarField magW ("magW", magSqr(Un-Us));

			gamma = 0.5*(fvc::ddt(rhon)+fvc::div(phin)-fvc::ddt(rhos)-fvc::div(phis));
			continuity = (fvc::ddt(rhon)+fvc::div(phin) + fvc::ddt(rhos)+ fvc::div(phis));

			#include "TEqn.H"
			#include "UEqns.H"

            // --- PISO loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

        }

        runTime.write();


		if(runTime.outputTime())
		{
			volVectorField GM ("GM", A*rhon*rhos*magSqr(Un-Us)*(Un-Us));
			GM.write();
		}
		
        
		Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
