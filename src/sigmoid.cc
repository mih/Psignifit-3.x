/*
 *   See COPYING file distributed along with the psignifit package for
 *   the copyright and license terms
 */
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

/************************************************************
 * Gumbel_l cdf
 */

double PsiGumbelL::f ( double x )
{
	if (x!=lastx) {
		lastx = x;
		lastf = 1-exp(-exp(x));
	}
	return lastf;
}

double PsiGumbelL::df ( double x )
{
	if ( x!=lastdx ) {
		lastdx = x;
		lastdf = exp( x - exp(x));
	}
	return lastdf;
}

double PsiGumbelL::ddf ( double x )
{
	if ( x!=lastddx ) {
		lastddx = x;
		lastddf = exp ( x - exp(x) ) * (1-exp(x));
	}
	return lastddf;
}

double PsiGumbelL::inv ( double p )
{
	if ( p!=lastp ) {
		lastp = p;
		lastinvp = log(-log(1-p));
	}
	return lastinvp;
}

/************************************************************
 * Gumbel_r cdf
 */

double PsiGumbelR::f ( double x )
{
	if (x!=lastx) {
		lastx = x;
		lastf = exp(-exp(-x));
	}
	return lastf;
}

double PsiGumbelR::df ( double x )
{
	if ( x!=lastdx ) {
		lastdx = x;
		lastdf = exp(-x-exp(-x));
	}
	return lastdf;
}

double PsiGumbelR::ddf ( double x )
{
	if ( x!=lastddx ) {
		lastddx = x;
		lastddf = exp ( -x - exp(-x) ) * (exp(-x)-1);
	}
	return lastddf;
}

double PsiGumbelR::inv ( double p )
{
	if ( p!=lastp ) {
		lastp = p;
		lastinvp = log(-log(p));
	}
	return lastinvp;
}

/************************************************************
 * Cauchy cdf
 */

double PsiCauchy::f ( double x )
{
	return atan ( x )/M_PI + 0.5;
}

double PsiCauchy::df ( double x )
{
	return 1./(M_PI*(1+x*x));
}

double PsiCauchy::ddf ( double x )
{
	return -2*x/( M_PI * (1+2*x*x+4*x*x*x*x) );
}

double PsiCauchy::inv ( double p )
{
	return tan ( M_PI*(p-0.5) );
}

/************************************************************
 * Exponential cdf
 */

double PsiExponential::f ( double x )
{
	if (x<0)
		return 0;
	else
		return 1-exp ( -x );
}

double PsiExponential::df ( double x )
{
	if (x<0)
		return 0;
	else
		return exp( -x );
}

double PsiExponential::ddf ( double x )
{
	if (x<0)
		return 0;
	else
		return -exp( -x );
}

double PsiExponential::inv ( double p )
{
	if ( p>0 && p<1 )
		return -log(1-p);
	else
		throw BadArgumentError();
}
