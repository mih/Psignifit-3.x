#ifndef PRIOR_H
#define PRIOR_H

#include <cmath>
#include "special.h"

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

/** \brief Uniform prior on an interval
 *
 * This prior defines a uniform distribution on an interval. It's pdf is thus \f$u-l\f$ if \f$x\in(l,u)\f$ and
 * 0 otherwise.
 */
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
 *
 * This defines a gaussian prior on the entire real axis. It's pdf is defined by the normal
 *
 \f[
 f(x) = \frac{1}{\sqrt{2 \pi}\sigma} \exp \left( - \frac{(x-\mu)^2}{2\sigma^2} \right).
 \f]
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
		GaussPrior ( double mean, double sd ) : mu(mean), sg(sd), normalization(1./(sqrt(2*M_PI)*sigma)), twovar(2*sg*sg), var(sg*sg) {}        ///< initialize prior to have mean mean and standard deviation sd
		double pdf ( double x ) { return normalization * exp ( - (x-mu)*(x-mu)/twovar ); }                                              ///< return pdf of the prior at position x
		double dpdf ( double x ) { return - x * pdf ( x ) / var; }                                                                      ///< return derivative of the prior at position x
};

/** \brief beta prior
 *
 * This defines a beta prior distribution that is defined on an interval. It's pdf is defined by
 *
 \f[
 f(x) = \frac{x^{\alpha-1} (1-x)^{\beta-1}} {B(\alpha,\beta)},
 \f]
 *
 * if \f$x\in(0,1)\f$. It is zero otherwise.
 */
class BetaPrior : public PsiPrior
{
	private:
		double alpha;
		double beta;
		double normalization;
	public:
		BetaPrior ( double al, double bt ) : alpha(al), beta(bt), normalization(betaf(al,bt)) {}                      ///< Initialize with parameters alpha=al, beta=bt
		double pdf ( double x ) { return (x<0||x>1 ? 0 : pow(x,alpha-1)*pow(1-x,beta-1)/normalization); }             ///< return beta pdf
		double dpdf ( double x ) { return (x<0||x>1 ? 0 : ((alpha-1)*pow(x,alpha-2) + (beta-1)*pow(1-x,beta-2))/normalization); }      ///< return derivative of beta pdf
};

/** \brief gamma prior
 *
 * This defines a gamma prior that is defined for the positive axis. It's pdf is defined by
 *
 \f[
 f(x) = x^{k-1} \frac{\exp(-x/\theta)}{\theta^k\Gamma(k)},
 \f]
 * for positive numbers and it is zero otherwise.
 */
class GammaPrior : public PsiPrior
{
	private:
		double k;
		double theta;
		double normalization;
	public:
		GammaPrior ( double shape, double scale ) : k(shape), theta(scale), normalization(pow(scale,shape)*exp(gammaln(shape))) {}                         ///< Initialize a gamma prior
		double pdf ( double x ) { return (x>0 ? pow(x,k-1)*exp(-x/theta)/normalization : 0 );}                                                             ///< return pdf at position x
		double dpdf ( double x ) { return (x>0 ? ( (k-1)*pow(x,k-2)*exp(-x/theta)-pow(x,k-1)*exp(-x/theta)/theta)/normalization : 0 ); }                   ///< return derivative of pdf
};

#endif
