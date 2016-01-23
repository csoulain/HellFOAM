fvVectorMatrix UsEqn(Us, Us.dimensions()*dimVol/dimTime*dimDensity);
fvVectorMatrix UnEqn(Un, Un.dimensions()*dimVol/dimTime*dimDensity);

{
    {
        UsEqn =
        (

            fvm::ddt(rhos,Us)
          + fvm::div(phis, Us, "div(phis,Us)")

          - fvm::laplacian(mus, Us)

		  - s*rhos*fvc::grad(T)


		  - aGM*A*rhos*rhon*magW*Un
		  + fvm::Sp(aGM*A*rhos*rhon*magW,Us)

		  + fvm::Sp(0.5*gamma,Us)
		  + 0.5*gamma*Us

		  - aGW*0.5*rhos*rhon/rho*fvc::grad(magW)

        );

        UsEqn.relax();
  //      UsEqn.solve();

    }

    {
        UnEqn =
        (
            fvm::ddt(rhon,Un)
          + fvm::div(phin, Un, "div(phin,Un)")

          - fvm::laplacian(mun, Un)
		  + rhos*s*fvc::grad(T)

		  - fvm::Sp(0.5*gamma,Un)
		  - 0.5*gamma*Us

		  - aGM*A*rhon*rhos*magW*Us
		  + fvm::Sp(aGM*A*rhon*rhos*magW,Un)

		  + aGW*0.5*rhon*rhos/rho*fvc::grad(magW)
	
        );

        UnEqn.relax();
   //     UnEqn.solve();

    }
}
