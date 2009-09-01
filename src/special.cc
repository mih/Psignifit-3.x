#include "special.h"

const double SQRT2PI ( sqrt(2*M_PI) );

double Phi ( double x ) {
	return 0.5*(1+erf(x/M_SQRT2));
}

double invPhi ( double p ) {
	/*
	double p0,p1,pp;
	double a,b;
	double x0,x1,x(0);
	double stepwidth;

	// There is no closed form solution
	// use secand method to solve
	do {
		stepwidth = 0.02*drand48();
		x0 = x-stepwidth; x1 = x+stepwidth;
		p0 = Phi(x0)-p; p1 = Phi(x1)-p;
		b = (p1-p0)/(x1-x0); a = p0-b*x0;
		x = -a/b;
	} while ( fabs(Phi(x)-p)>1e-7);

	return x;
	// return M_SQRT2*erfinv(2*p-1);
	*/
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
