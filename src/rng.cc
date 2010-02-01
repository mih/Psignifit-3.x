/*
 *   See COPYING file distributed along with the psignifit package for
 *   the copyright and license terms
 */
#include "rng.h"

double GaussRandom::draw ( void )
{
	if ( good ) {
		good = 0;
		return y*sigma+mu;
	} else {
		do {
			x1 = 2*rngcall () - 1;
			x2 = 2*rngcall () - 1;
			w = x1*x1 + x2*x2;
		} while ( w >= 1.0 );

		w = sqrt ( -2.*log(w)/w );
		y = x2*w;
		good = true;
		return x1*w*sigma+mu;
	}
}

double BinomialRandom::draw ( void )
{
	int k(0),i;
	for (i=0; i<n; i++)
		if ( rngcall() < p )
			k++;
	return k;
}
