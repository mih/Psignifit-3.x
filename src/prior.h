#ifndef PRIOR_H
#define PRIOR_H

#include <cmath>

/** \brief base class for all priors
 *
 * This default prior does nothing in particular. It poses no restriction on the respective parameter at all.
 * at any value x, the prior returns 1. Thus it is an "improper" prior in the sense that it does not correspond
 * to a proper probability distribution.
 */
class PsiPrior
{
	public:
		virtual double pdf ( double x ) { return 1.;}    ///< evaluate the pdf of the prior at position x (in this default form, the parameter is completely unconstrained)
		virtual double dpdf ( double x ) { return 0.; }  ///< evaluate the derivative of the pdf of the prior at position x (in this default form, the parameter is completely unconstrained)
};

/** \brief Uniform prior on an interval */
class UniformPrior : public PsiPrior
{
	private:
		double lower;
		double upper;
		double height;
	public:
		UniformPrior ( double low, double high ) : lower(low), upper(high), height(1./(high-low)) {} ///< Set up a UniformPrior on the interval from low to high
		double pdf ( double x ) { return ( x>lower && x<upper ? height : 0 ); }                      ///< evaluate the pdf of the prior at position x
		double dpdf ( double x ) { return ( x!=lower && x!=upper ? 0 : (x==lower ? 1e20 : -1e20 ));} ///< derivative of the pdf of the prior at position x (jumps at lower and upper are replaced by large numbers)
};

/** \brief gaussian (normal) prior
 */
class GaussPrior : public PsiPrior
{
	private:
		double mu;
		double sg;
		double normalization;
		double var;
		double twovar;
	public:
		GaussPrior ( double mean, double sd ) : mu(mean), sg(sd), normalization(1./sqrt(2*M_PI)), twovar(2*sg*sg), var(sg*sg) {}        ///< initialize prior to have mean mean and standard deviation sd
		double pdf ( double x ) { return normalization * exp ( - (x-mu)*(x-mu)/twovar ); }                                              ///< return pdf of the prior at position x
		double dpdf ( double x ) { return - x * pdf ( x ) / var; }                                                                      ///< return derivative of the prior at position x
};

#endif
