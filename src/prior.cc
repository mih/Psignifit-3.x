/*
 *   See COPYING file distributed along with the psignifit package for
 *   the copyright and license terms
 */
#include "prior.h"

#include <iostream>

double BetaPrior::rand ( void )
{
	double x ( rng.draw() );
	double p ( rng.draw() );

	while ( p > (pdf(x)/mode) ) {
		x = rng.draw();
		p = rng.draw();
	}

	return x;
}

double GammaPrior::rand ( void )
{
	double al ( k-1 );
	double v ( -log(rng.draw()) );
	double w ( -log(rng.draw()) );

	while ( w < al*(v-log(v)-1) ) {
		v = -log( rng.draw() );
		w = -log( rng.draw() );
	}

	return k*v*theta;
}
