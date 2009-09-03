#include "special.h"

const double SQRT2PI ( sqrt(2*M_PI) );

double Phi ( double x ) {
	return 0.5*(1+erf(x/M_SQRT2));
}

double invPhi ( double p ) {
	double x(0),step;
	double g,gprime;

	do {
		g = Phi(x)-p;
		gprime = exp(-0.5*x*x)/SQRT2PI;
		step = g/gprime;
		x -= step;
	} while (fabs(step)>1e-7);

	return x;
}


double safe_log ( double x )
{
	return (x>0 ? log(x) : -1e20);
}

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

double betaf(double z, double w) {
	return exp(gammaln(z)+gammaln(w)-gammaln(z+w));
}
