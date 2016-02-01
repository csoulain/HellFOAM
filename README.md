# HellFOAM

Description

    Solve the He(II) superfluid equations. The solver is based on the Super-PISO
	algorithm. 

    If you use this program, please cite

	Soulaine, C.; Quintard, M.; Allain, H.; Baudouy, B. & Van Weelderen, R. 
	A PISO-like algorithm to simulate superfluid helium flow with the two-fluid model 
	Computer Physics Communications , 2015, 187, 20-28

	Soulaine, C.; Quintard, M.; Allain, H.; Baudouy, B. & Van Weelderen, R. 
	Numerical Investigation of Heat Transfer In a Forced Flow of He II 
	Proceedings of the 15th International Heat Transfer Conference, 
	IHTC-15 August 10-15, 2014, Kyoto, Japan, 2014

Tutorials
	
	The package comes with different examples (just run HellFoam in a terminal, the grid and 
	initialization are already set up)
		* capillaryGorterMellink: solve superfluid helium flow in a tube in the Gorter-Mellink
			regime. Simulation results can be compared with analytical solutions
			(see Soulaine et al. 2015, section 4.2)
		* capillaryWithCylindersGorterMellink: solve superfluid helium flow in a porous medium 
			model (a capillary containing several solid cylinders) in the Gorter-Mellink regime
			(see Soulaine et al. 2015, section 4.4)

Authors

	2013-07-25: Cyprien Soulaine (CS, cyprien.soulaine@gmail.com) - First version
	2013-09-09: CS - Second major release
					* implementation of temperature-dependant variables. 
					* Correction of the temperature equation.
					* Upgrade to OpenFOAM 2.2.1
	2014-03-24: CS - minor release
					* Upgrade to OpenFOAM 2.3.0 
	2016-01-08: CS - minor release
					* Upgrade to OpenFOAM 3.0.1 
