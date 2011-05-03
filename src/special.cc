/*
 *   See COPYING file distributed along with the psignifit package for
 *   the copyright and license terms
 */
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

double gammainc ( double x, double a ) {
	double gm,Ga,xi(1),h,dl,aa,bb;
	unsigned int i;
	double D,C;
	gm = exp ( -x ) * pow(x,a);
	Ga = exp ( gammaln(a) );
	if ( x < a+1 ) {
		// Series
		i = 0;
		while (i<2000) {
			gm += xi * Ga / ( exp ( gammaln ( a+1+i ) ) );
			xi *= x;
			i++;
		}
		return gm;
	} else {
		// Continued fraction
		h = 1e-30;
		D = x+1-a;
		C = x+1-a + 1/h;
		dl = C/D;
		h *= dl;
		i = 1;
		while ( i++<2000 ) {
			bb = x+2*(i-1)+1-a;
			aa = (i-1)*(i-1-a);
			D = bb + aa * D;
			C = bb + aa / C;
			if (D==0) D = 1e-30;
			if (C==0) C = 1e-30;
			dl = C/D;
			h *= dl;
			if ( fabs(dl-1) < 1e-7 ) break;
		}
		gm *= h;
		gm = Ga - gm;
	}
	return gm;
}

double betaf(double z, double w) {
	return exp(gammaln(z)+gammaln(w)-gammaln(z+w));
}

double betainc ( double x, double a, double b ) {
	unsigned int i,m;
	double h(1), C(1), D(0), d, dl;

	i = 0;
	while (i++<2000) {
		if ( i%2==0 ) {
			m = i/2;
			d = m*(b-m)*x / ((a+2*m-1)*(a+2*m));
		} else {
			m = (i-1)/2;
			d = - (a+m)*(a+b+m)*x / ( (a+2*m-1)*(a+2*m) );
		}

		D = 1 + d*D;
		C = 1 + d/C;
		if ( D==0 ) d = 1e-30;
		if ( C==0 ) d = 1e-30;

		dl = C/D;
		h = dl*h;
		if ( fabs(dl-1) < 1e-7 ) break;
	}

	return pow ( x,a )*pow ( x, b ) / (a*betaf ( a, b ) ) * h;
}

double psi ( double z ) {
	/* This algorithm is based on two identities:
	 * 1. The first is an approximation formula that works for large z
	 *       psi(z) ~ log(z) - 1./(2*z) - 1./12*z**2 + 1./120*z**4 - 1./252*z**6 + O(1./z**8)
	 *	See for example Abramowitz & Stegun eq 6.3.18
	 * 2. The second is a recursion formula (Abramowitz & Stegun eq 6.3.5)
	 *	     psi(z+1) = psi(z) + 1./z
	 * Thus, for sufficiently large z (z>6 results in accuracies of order 1e-9) we take the approximation
	 * formula. For smaller z, we apply the second formula as a telescope.
	 */
	if ( z > 6 ) {
		return log ( z ) - 1./(2*z) - 1./(12*z*z) + 1./(120*z*z*z*z) - 1./(252*z*z*z*z*z*z);
	} else {
		return psi ( z+1 ) - 1./z;
	}
}

double digamma ( double z ) {
	/*
	 * This algorithm is similar to that for evaluation of the psi function.
	 * 1. There is again an approximation formula for large z
	 *		digamma(z) ~ 1./z + 1./(2*z**2) - 1./(6*z**3) + 1./(42*z**5) - 1./(30*z**9) + O(1./z**10)
	 *	See for example Abramowitz & Stegun eq 6.4.12
	 * 2. There is a recursion formula (Abraomowitz & Stegun eq 6.4.6
	 *      digamma(z) ~ digamma(z+1) = digamma(z) - 1./z**2
	 * Again, we combine these two and use the recursion for z<5. This results in accuracies of order 1e-9.
	 */
	if ( z > 5 )
		return 1./z + 1./(2*z*z) + 1./(6*z*z*z) - 1./(30*z*z*z*z*z) + 1./(42*z*z*z*z*z*z*z) - 1./(30*z*z*z*z*z*z*z*z*z);
	else
		return digamma ( z+1 ) + 1./(z*z);
}
