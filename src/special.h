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

#endif
