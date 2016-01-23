
//surfaceScalarField phiT ("phiT", fvc::interpolate(rho*s/T)*phin);

surfaceScalarField phinOverRhon ("phinT", phin/fvc::interpolate(rhon));


fvScalarMatrix TEqn
(
	(s*dRhodT + rho*dSdT)*(fvm::ddt(T) + fvm::div(phinOverRhon,T,"div(phi,T)") ) - fvm::laplacian(kf/T,T, "laplacian(DT,T)") 
	== 
	  Q/fvc::domainIntegrate(alpha)*alpha/T  + aGM*A*rhon*rhos*magW*magW/T
	- fvm::Sp((rho*s/T - (s*dRhodT + rho*dSdT))*fvc::div(phinOverRhon),T)
);


TEqn.relax();
TEqn.solve();
