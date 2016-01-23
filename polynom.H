

dimensionedScalar calcPoly( List<scalar> const & polyCoeff, dimensionedScalar const T, dimensionedScalar & polyT)
{
// calcule P(T)=a0 + a1*T + a2*T^2 + ....... + an*T^n 
// où T est un champ scalair et (ai)i est le polynome

	polyT.value() = 0.0;			

	for(int i=0; i<polyCoeff.size(); i++)
	{
		polyT.value() += polyCoeff[i]*Foam::pow(T.value(),i);			
	}

	return polyT;
}


void printPoly(  List<scalar> const & polyCoeff )
{
	for(int k=0; k<polyCoeff.size() ; k++)
	{
		Info << "Poly["<<k<<"] = " << polyCoeff[k] << endl;
	}		
	Info << nl << endl;
}


List<scalar> derivPoly1( List<scalar> const & polyCoeff)
{
	List<scalar> derivPolyCoeff(polyCoeff.size()-1);

	derivPolyCoeff = 0;

	for(int k=0; k<derivPolyCoeff.size() ; k++)
	{
		derivPolyCoeff[k] = (k+1)*polyCoeff[k+1];
//		Info << "D poly D T = " << derivPolyCoeff[k] << endl;
	}

	return derivPolyCoeff;
}


List<scalar> addPoly(List<scalar> const & polyCoeff1, List<scalar> const & polyCoeff2)
{
	// attention, lorsque les degrés des Polynomes sont différents, l'algo suivant
	// est peut être source d'erreurs.
	// attention aussi lors que les termes de plus haut degré s'annulent.

	List<scalar> polyCoeffResult(max(polyCoeff1.size(),polyCoeff2.size()));

	for(int k=0; k<polyCoeffResult.size() ; k++)
	{
		polyCoeffResult[k]=polyCoeff1[k]+polyCoeff2[k];
	}

	return polyCoeffResult;
}


List<scalar> multiplyPoly(List<scalar> const & polyCoeff1, List<scalar> const & polyCoeff2)
{
	// Multiplie 2 polynomes. Le résultat est la liste des coefficients du
	// nouveau polynome.

	List<scalar> polyCoeffResult(polyCoeff1.size()+polyCoeff2.size()-1);
	List<scalar> polyCoeff1_ = polyCoeff1;
	List<scalar> polyCoeff2_ = polyCoeff2;

	for(int k=polyCoeff1_.size() ; k<polyCoeffResult.size() ; k++)
	{
		polyCoeff1_.append(0);
	}

	for(int k=polyCoeff2_.size() ; k<polyCoeffResult.size() ; k++)
	{
		polyCoeff2_.append(0);
	}

	for(int k=0; k<polyCoeffResult.size() ; k++)
	{
		polyCoeffResult[k]=0.0;
		for(int m=0 ; m<=k ; m++)
		{
			polyCoeffResult[k]+=polyCoeff1_[m]*polyCoeff2_[k-m];
		}
		
	}
	polyCoeff1_.clear();
	polyCoeff2_.clear();

	return polyCoeffResult;
}


