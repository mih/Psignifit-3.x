#ifndef SIGMOID_H
#define SIGMOID_H

#include "errors.h"
#include "special.h"
#include <cmath>

/** \brief common base class for all sigmoids */
class PsiSigmoid
{
	public:
		virtual double f   ( double x ) { throw NotImplementedError(); }      ///< This should return the value of the sigmoid itself (between 0 and 1)
		virtual double df  ( double x ) { throw NotImplementedError(); }      ///< This should give the first derivative of the sigmoid
		virtual double ddf ( double x ) { throw NotImplementedError(); }      ///< This should give the second derivative of the sigmoid
		virtual double inv ( double p ) { throw NotImplementedError(); }      ///< This should give the inverse of the sigmoid (taking values between 0 and 1)
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
};

#endif
