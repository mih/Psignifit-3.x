/*
 *   See COPYING file distributed along with the psignifit package for
 *   the copyright and license terms
 */
#include "rng.h"

long __psignifit_rng_idum (-1);

double PsiRandom::rngcall ( void ) {
	// Long period (>2*10e18) random number generator of L'Ecuyer with Bays-Durham shuffle
	// and added safeguards. Returns a uniform random deviate between 0.0 and 1.0 (exclusive
	// of the endpoints). Uses __psignifit_rng_idum as seed which should be negative and not
	// be changed during a simulation run.
	// RNMX should approximate the largest floating point value that is less than 1.
	//
	// This has been copied from numerical recipes
	long *idum = &__psignifit_rng_idum;
	int j;
	long k;
	static long idum2 ( 123456789 );
	static long iy ( 0 );
	static long iv[NTAB];
	double temp;

	if ( *idum <= 0 ) {
		if ( -(*idum) < 1 ) *idum=1;
		else *idum = (-*idum);
		idum2 = (*idum);
		for ( j = NTAB+7; j>=0; j-- ) {
			k = (*idum) / IQ1;
			*idum = IA1*(*idum-k*IQ1) - k*IR1;
			if ( *idum<0 ) *idum += IM1;
			if ( j<NTAB ) iv[j] = *idum;
		}
		iy = iv[0];
	}
	k = (*idum)/IQ1;
	*idum = IA1*(*idum-k*IQ1) - k*IR1;
	if ( *idum<0 ) *idum += IM1;
	k = idum2 / IQ2;
	idum2 = IA2*(idum2-k*IQ2) - k*IR2;
	j = iy/NDIV;
	iy = iv[j] - idum2;
	iv[j] = *idum;
	if ( iy<1 ) iy += IMM1;
	if ( (temp=AM*iy)>RNMX ) return RNMX;
	else return temp;
}

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

void setSeed(long int seedval){
	if ( seedval==0 ) {
		__psignifit_rng_idum = -1;
	} else if ( seedval > 0 ) {
		__psignifit_rng_idum = -seedval;
	} else {
		__psignifit_rng_idum = seedval;
	}
}
