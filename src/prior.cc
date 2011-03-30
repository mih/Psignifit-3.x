/*
 *   See COPYING file distributed along with the psignifit package for
 *   the copyright and license terms
 */
#include "prior.h"

#include <iostream>

void GaussPrior::shrink ( double xmin, double xmax ) {
	double s ( 0.5*(xmax-xmin) ), m ( 0.5*(xmin+xmax) );
	if ( s<std() ) {
		mu = m;
		sg = s;
		var = sg*sg;
		twovar = 2*var;
		rng = GaussRandom ( mu, sg );
		normalization = 1./(sqrt(2*M_PI)*sg);
	}
}

void BetaPrior::shrink ( double xmin, double xmax ) {
	double s ( 0.5*(xmax-xmin) ), m ( 0.5*(xmin+xmax) );
	if ( s<std() ) {
		beta = m*(1-m)*(1-m)/(s*s) - 1 + m;
		alpha = m*beta/(1-m);
		normalization = betaf(alpha,beta);
		rng = BetaRandom ( alpha, beta );
	}
}

void GammaPrior::shrink ( double xmin, double xmax ) {
	double xr ( xmin/xmax );
	k = ((1+xr)/(1-xr));
	k *= k;
	theta = xmax / (k+sqrt(k));
	normalization = pow(theta,k)*exp(gammaln(k));
	rng = GammaRandom ( k, theta );
}

void nGammaPrior::shrink ( double xmin, double xmax ) {
	double ymin(-xmax), ymax(-xmin);
	GammaPrior::shrink ( ymin, ymax );
}
