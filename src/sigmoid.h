/*
 *   See COPYING file distributed along with the psignifit package for
 *   the copyright and license terms
 */
#ifndef SIGMOID_H
#define SIGMOID_H

#include "errors.h"
#include "special.h"
#include <cmath>

/** \brief common base class for all sigmoids */
class PsiSigmoid
{
	public:
		virtual double f   ( double x ) { throw NotImplementedError(); }            ///< This should return the value of the sigmoid itself (between 0 and 1)
		virtual double df  ( double x ) { throw NotImplementedError(); }            ///< This should give the first derivative of the sigmoid
		virtual double ddf ( double x ) { throw NotImplementedError(); }            ///< This should give the second derivative of the sigmoid
		virtual double inv ( double p ) { throw NotImplementedError(); }            ///< This should give the inverse of the sigmoid (taking values between 0 and 1)
		virtual int    getcode ( void ) const { throw NotImplementedError(); }            ///< return the sigmoid identifier
};

/** \brief logistic function
 *
 * The logistic function is given by f(x) = 1/(1+exp(-x))
 */
class PsiLogistic : public PsiSigmoid
{
	private:
		double lastx;
		double lastfx;
	public:
		PsiLogistic ( void ) : lastx(1e20) {}  ///< constructor
		double f ( double x );                 ///< value of the sigmoid at position x
		double df ( double x );                ///< derivative of the sigmoid at position x
		double ddf ( double x );               ///< second derivative of the sigmoid
		double inv ( double p ) { return log(p/(1-p)); }  ///< inverse of the sigmoid
		int getcode ( void ) const { return 1; }     ///< return the sigmoid identifier
};

/** \brief gaussian cdf function
 *
 * The gaussian cdf function is given by f(x) = Phi(x), where Phi is the cumulative distribution function for the gaussian distribution.
 */
class PsiGauss : public PsiSigmoid
{
	private:
		double lastx;
		double lastf;
		double lastx_d;
		double lastdf;
		double lastx_dd;
		double lastddf;
		double lastp;
		double lastinvp;
	public:
		PsiGauss ( void ) : lastx(1e20), lastx_d(1e20),lastx_dd(1e20),lastp(1e20) {};
		double f   ( double x );                 ///< value of the sigmoid at x
		double df  ( double x );                 ///< derivative of the sigmoid at x
		double ddf ( double x );                 ///< second derivative of the sigmoid at x
		double inv ( double p );                 ///< inverse of the sigmoid
		int getcode ( void ) const { return 2; }       ///< return the sigmoid identifier
};

/** \brief left-skewed gumbel cdf
 *
 * Cumulative densitiy function of the gumbel distribution.
 */
class PsiGumbelL : public PsiSigmoid
{
	private:
		double lastx;
		double lastf;
		double lastdx;
		double lastddx;
		double lastdf;
		double lastddf;
		double lastp;
		double lastinvp;
	public:
		PsiGumbelL ( void ) : lastx(-1e20), lastdx(-1e20), lastddx(-1e20), lastp(0), lastinvp(-1e20) {}
		double f   ( double x );              ///< returns the value of the gumbel cdf at position x
		double df  ( double x );              ///< returns the derivative of the gumbel cdf at position x
		double ddf ( double x );              ///< returns the 2nd derivative of the gumbel cdf at position x
		double inv ( double p );              ///< returns the inverse of the gumbel cdf at position p
		int getcode ( void ) const { return 3; }    ///< return the sigmoid identifier
};

/** \brief right-skewed gumbel cdf
 *
 * Cumulative densitiy function of the gumbel distribution.
 */
class PsiGumbelR : public PsiSigmoid
{
	private:
		double lastx;
		double lastf;
		double lastdx;
		double lastddx;
		double lastdf;
		double lastddf;
		double lastp;
		double lastinvp;
	public:
		PsiGumbelR ( void ) : lastx(-1e20), lastdx(-1e20), lastddx(-1e20), lastp(0), lastinvp(-1e20) {}
		double f   ( double x );             ///< returns the value of the right skewed gumbel cdf at position x
		double df  ( double x );             ///< returns the derivative of the right skewed gumbel cdf at position x
		double ddf ( double x );             ///< returns the 2nd derivative of the right skewed gumbel cdf at position x
		double inv ( double p );             ///< returns the inverse of the right skewed gumbel cdf at position p
		int getcode ( void ) const { return 3; }   ///< return the sigmoid identifier
};

/** \brief cauchy cdf
 *
 * Cumulative density function of the cauchy distribution
 */
class PsiCauchy : public PsiSigmoid
{
	public:
		double f   ( double x );             ///< returns the value of the cauchy cdf at position x
		double df  ( double x );             ///< returns the derivative of the cauchy cdf at position x
		double ddf ( double x );             ///< returns the 2nd derivative of the cauchy cdf at position x
		double inv ( double p );             ///< returns the inverse of the cauchy cdf at position x
		int    getcode ( void ) const { return 4; }///< returns the sigmoid identifier
};

/** \brief exponential cdf
 *
 * Cumulative density function of the exponential distribution
 * combined with a polyCore this will give a weibull
 */
class PsiExponential : public PsiSigmoid
{
	public:
		double f   (double x );              ///< returns the value of the exponential cdf at position x
		double df  (double x );              ///< returns the derivative of the exponential cdf at position x
		double ddf (double x );              ///< returns the 2nd derivative of the exponential cdf at position x
		double inv (double p );              ///< returns the return the inverse of the exponential cdf at position x
		int    getcode ( void ) const { return 5; }///< returns the sigmoid identifier
};

#endif
