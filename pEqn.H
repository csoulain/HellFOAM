{
    volScalarField rUsA(1.0/UsEqn.A());
    volScalarField rUnA(1.0/UnEqn.A());

    phis == fvc::interpolate(rhos)*(fvc::interpolate(Us) & mesh.Sf());
    phin == fvc::interpolate(rhon)*(fvc::interpolate(Un) & mesh.Sf());

    surfaceScalarField rUsAf(fvc::interpolate(rhos*rUsA));
    surfaceScalarField rUnAf(fvc::interpolate(rhon*rUnA));

    Us = rUsA*UsEqn.H();
    Un = rUnA*UnEqn.H();


    phis = ((fvc::interpolate(rhos*Us) & mesh.Sf()) + rUsAf*fvc::ddtCorr(rhos, Us, phis));

    phin = ((fvc::interpolate(rhon*Un) & mesh.Sf()) + rUnAf*fvc::ddtCorr(rhon, Un, phin));

//	surfaceScalarField phiTot ("phiTot", fvc::interpolate(rhos)*phis+fvc::interpolate(rhon)*phin);

	surfaceScalarField phiTot ("phiTot", phis+phin);

    surfaceScalarField Dp
    (
        "Dp",
    	    fvc::interpolate(rhos)*rUsAf/fvc::interpolate(rho)
		 +  fvc::interpolate(rhon)*rUnAf/fvc::interpolate(rho)
    );

    for (int nonOrth=0; nonOrth<=pimple.nNonOrthCorr(); nonOrth++)
    {
        fvScalarMatrix pEqn
        (
            fvm::laplacian(Dp, p) == fvc::div(phiTot) + fvc::ddt(rho)
        );

        pEqn.setReference(pRefCell, pRefValue);

/*
        pEqn.solve
        (
            mesh.solver(p.select(pimple.finalInnerIter(corr, nonOrth)))
        );
*/
		pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));

        if (nonOrth == pimple.nNonOrthCorr())
        {
            surfaceScalarField SfGradp(pEqn.flux()/Dp);

            phis -= rUsAf*SfGradp*fvc::interpolate(rhos)/fvc::interpolate(rho);

            phin -= rUnAf*SfGradp*fvc::interpolate(rhon)/fvc::interpolate(rho);

			phi =  phin + phis;
       		//phi =  fvc::interpolate(rhon/rho)*phin + fvc::interpolate(rhos/rho)*phis;

            p.relax();
            SfGradp = pEqn.flux()/Dp;

//            Us += fvc::reconstruct(- rUsAf*SfGradp/fvc::interpolate(rho));
            Us += fvc::reconstruct(- rUsAf*SfGradp/fvc::interpolate(rho));
            Us.correctBoundaryConditions();

//            Un += fvc::reconstruct(- rUnAf*SfGradp/fvc::interpolate(rho));
            Un += fvc::reconstruct(- rUnAf*SfGradp/fvc::interpolate(rho));
            Un.correctBoundaryConditions();

            U = (rhos*Us + rhon*Un)/rho;
        }
    }
}

#include "continuityErrs.H"
