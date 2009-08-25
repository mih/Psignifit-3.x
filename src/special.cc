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
