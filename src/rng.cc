#include "rng.h"

double GaussRandom::draw ( void )
{
	if ( good ) {
		good = 0;
		return y*sigma+mu;
	} else {
		do {
			x1 = 2*drand48 () - 1;
			x2 = 2*drand48 () - 1;
			w = x1*x1 + x2*x2;
		} while ( w >= 1.0 );

		w = sqrt ( -2.*log(w)/w );
		y = x2*w;
		good = true;
		return x1*w*sigma+mu;
	}
}
