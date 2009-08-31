#ifndef SPECIAL_H
#define SPECIAL_H

#include <cmath>
#include <cstdlib>

/** \brief gaussian cumulative distribution function */
double Phi ( double x );

/** \brief inverse of the gaussian cumulative distribution function
 *
 * This function is really expensive. It inverts the gaussian cdf by performing a numerical
 * solution of the equation
 * Phi(x) - p = 0
 */
double invPhi ( double p );

/** \brief logarithm that does not return nan but instead a very low value (-1e20) */
double safe_log ( double x );

/** \brief logarithm of the gamma function */
double gammaln(double xx) {
	// More or less copied from Numerical Recipes
	double x,y,tmp,ser;
	static double cof[6]={
		76.18009172947146,
		-86.50532032941677,
		24.01409824083091,
		-1.231739572450155,
		0.1208650973866179e-2,
		-0.5395239384953e-5 };
	int j;

	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0; j<=5; j++) ser += cof[j]/++y;

	return -tmp+log(2.5066282746310005*ser/x);
}

/** beta function */
double betaf(double z, double w) {
	return exp(gammaln(z)+gammaln(w)-gammaln(z+w));
}

#endif
