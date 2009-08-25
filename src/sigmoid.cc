#include "sigmoid.h"

#ifdef DEBUG_SIGMOID
#include <iostream>
#endif

/************************************************************
 * Logistic sigmoid
 */

double PsiLogistic::f ( double x )
{
#ifdef DEBUG_SIGMOID
	std::cerr << "In logistic function.\n";
#endif
	if (x!=lastx) {
		lastx = x;
		lastfx = 1./(1.+exp(-x));
	}
#ifdef DEBUG_SIGMOID
	std::cerr << " lastx = " << lastx << "\n lastfx = " << lastfx << "\n";
#endif
	return lastfx;
}

double PsiLogistic::df ( double x )
{
	return f(x)*(1-f(x));
}

double PsiLogistic::ddf ( double x )
{
	return f(x)*(1-f(x))*(1-2*f(x));
}

/************************************************************
 * Gauss-CDF
 */

double PsiGauss::f ( double x )
{
	if (x==lastx)
		return lastf;
	else {
		lastx = x;
		lastf = Phi(x);
		return lastf;
	}
}

double PsiGauss::df ( double x )
{
	if (x==lastx_d)
		return lastdf;
	else {
		lastx_d = x;
		lastdf = exp ( - 0.5*x*x ) / sqrt(2*M_PI);
		return lastdf;
	}
}

double PsiGauss::ddf ( double x )
{
	if (x==lastx_dd)
		return lastddf;
	else {
		lastx_dd = x;
		lastddf = -x*df(x);
		return lastddf;
	}
}

double PsiGauss::inv ( double p )
{
	if (p==lastp)
		return lastinvp;
	else {
		lastp = p;
		lastinvp = invPhi(p);
		return lastinvp;
	}
}
