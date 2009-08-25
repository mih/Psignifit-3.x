#ifndef PRIOR_H
#define PRIOR_H

/** \brief base class for all priors
 *
 * This default prior does nothing in particular. It poses no restriction on the respective parameter at all.
 * at any value x, the prior returns 1. Thus it is an "improper" prior in the sense that it does not correspond
 * to a proper probability distribution.
 */
class PsiPrior
{
	public:
		virtual double pdf ( double x ) { return 1.;}  ///< evaluate the pdf of the prior at position x (in this default form, the parameter is completely unconstrained)
};

class UniformPrior : public PsiPrior
{
	private:
		double lower;
		double upper;
		double height;
	public:
		UniformPrior ( double low, double high ) : lower(low), upper(high), height(1./(high-low)) {} ///< Set up a UniformPrior on the interval from low to high
		double pdf ( double x ) { return ( x>lower && x<upper ? height : 0 ); }                      ///< evaluate the pdf of the prior at position x
};

#endif
