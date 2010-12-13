/*
 *   See COPYING file distributed along with the psignifit package for
 *   the copyright and license terms
 */
#ifndef PRIOR_H
#define PRIOR_H

#include <cmath>
#include "rng.h"
#include "special.h"

/** \brief base class for all priors
 *
 * This default prior does nothing in particular. It poses no restriction on the respective parameter at all.
 * at any value x, the prior returns 1. Thus it is an "improper" prior in the sense that it does not correspond
 * to a proper probability distribution.
 */
class PsiPrior
{
	private:
		PsiRandom rng;
	public:
		virtual double pdf ( double x ) { return 1.;}    ///< evaluate the pdf of the prior at position x (in this default form, the parameter is completely unconstrained)
		virtual double dpdf ( double x ) { return 0.; }  ///< evaluate the derivative of the pdf of the prior at position x (in this default form, the parameter is completely unconstrained)
		virtual double rand ( void ) { return rng.draw(); } ///< draw a random number
		virtual PsiPrior * clone ( void ) const { throw NotImplementedError(); }///< clone by value
		virtual double mean ( void ) const { return 0; } ///< return the mean
		virtual double std  ( void ) const { return 100; } ///< return the standard deviation
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
		UniformRandom rng;
	public:
		UniformPrior ( double low, double high ) : lower(low), upper(high), rng ( low, high ) { height = 1./(high-low); } ///< Set up a UniformPrior on the interval from low to high
		UniformPrior ( const UniformPrior& original ) : lower(original.lower),
                                                        upper(original.upper),
                                                        height(original.height),
                                                        rng(original.rng) {} ///< copy constructor
		double pdf ( double x ) { return ( x>lower && x<upper ? height : 0 ); }                      ///< evaluate the pdf of the prior at position x
		double dpdf ( double x ) { return ( x!=lower && x!=upper ? 0 : (x==lower ? 1e20 : -1e20 ));} ///< derivative of the pdf of the prior at position x (jumps at lower and upper are replaced by large numbers)
		double rand ( void ) { return rng.draw(); }                                                 ///< draw a random number
        PsiPrior * clone ( void ) const { return new UniformPrior(*this); }
		double mean ( void ) const { return 0.5*(lower+upper); }  ///< return the mean
		double std  ( void ) const { return sqrt((upper-lower)*(upper-lower)/12); }
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
		GaussRandom rng;
	public:
		GaussPrior ( double mean, double sd ) : mu(mean), sg(sd), var(sg*sg), twovar(2*sg*sg), rng(mean,sd) { normalization = 1./(sqrt(2*M_PI)*sg); }        ///< initialize prior to have mean mean and standard deviation sd
		GaussPrior ( const GaussPrior& original ) : mu(original.mu),
                                                    sg(original.sg),
                                                    normalization(original.normalization),
                                                    var(original.var),
                                                    twovar(original.twovar),
                                                    rng(original.rng) {} ///< copy contructor
		double pdf ( double x ) { return normalization * exp ( - (x-mu)*(x-mu)/twovar ); }                                              ///< return pdf of the prior at position x
		double dpdf ( double x ) { return - x * pdf ( x ) / var; }                                                                      ///< return derivative of the prior at position x
		double rand ( void ) {return rng.draw(); }
        PsiPrior * clone ( void ) const { return new GaussPrior(*this); }
		double mean ( void ) const { return mu; } ///< mean
		double std  ( void ) const { return sg; } ///< return standard deviation
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
		UniformRandom rng;
		double mode;
	public:
		BetaPrior ( double al, double bt ) : alpha(al), beta(bt), normalization(betaf(al,bt)), rng (0,1) {
            mode = (al-1)/(al+bt-2);
            mode = pdf(mode); }                      ///< Initialize with parameters alpha=al, beta=bt
        BetaPrior ( const BetaPrior& original) : alpha(original.alpha),
                                                 beta(original.beta),
                                                 normalization(original.normalization),
                                                 rng(original.rng),
                                                 mode(original.mode) {} ///< copy constructor
		double pdf ( double x ) { return (x<0||x>1 ? 0 : pow(x,alpha-1)*pow(1-x,beta-1)/normalization); }             ///< return beta pdf
		double dpdf ( double x ) { return (x<0||x>1 ? 0 : ((alpha-1)*pow(x,alpha-2)*pow(1-x,beta-1) + (beta-1)*pow(1-x,beta-2)*pow(x,alpha-1))/normalization); }      ///< return derivative of beta pdf
		double rand ( void );                                                                                         ///< draw a random number using rejection sampling
        PsiPrior * clone ( void ) const { return new BetaPrior(*this); }
		double mean ( void ) const { return alpha/(alpha+beta); }
		double std  ( void ) const { return sqrt ( alpha*beta/((alpha+beta)*(alpha+beta)*(alpha+beta+1)) ); }
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
		UniformRandom rng;
	public:
		GammaPrior ( double shape, double scale ) : k(shape), theta(scale) { normalization = pow(scale,shape)*exp(gammaln(shape)); }                         ///< Initialize a gamma prior
        GammaPrior ( const GammaPrior& original ) : k(original.k),
                                                    theta(original.theta),
                                                    normalization(original.normalization),
                                                    rng(original.rng) {} ///< copy constructor
		virtual double pdf ( double x ) { return (x>0 ? pow(x,k-1)*exp(-x/theta)/normalization : 0 );}                                                             ///< return pdf at position x
		virtual double dpdf ( double x ) { return (x>0 ? ( (k-1)*pow(x,k-2)*exp(-x/theta)-pow(x,k-1)*exp(-x/theta)/theta)/normalization : 0 ); }                   ///< return derivative of pdf
		virtual double rand ( void );
        PsiPrior * clone ( void ) const { return new GammaPrior(*this); }
		virtual double mean ( void ) const { return k*theta; }
		double std  ( void ) const { return sqrt ( k*theta*theta ); }
};

/** \brief negative gamma prior
 *
 * This defines a gamma prior that is defined for the negative axis. Thus, the pdf is defined by
 *
 \f[
 f(x) = (-x)^{k-1} \frac{\exp(x/\theta)}{\theta^k\Gamma(k)},
 \f]
 * for negative numbers and is zero otherwise.
 */
class nGammaPrior : public GammaPrior
{
	public:
		nGammaPrior ( double shape, double scale ) : GammaPrior(shape,scale) {}
        nGammaPrior ( const nGammaPrior& original ) : GammaPrior(original) {} ///< copy constructor
		double pdf ( double x ) { return GammaPrior::pdf ( -x ); }
		double dpdf ( double x ) { return GammaPrior::dpdf ( -x ); }
		double rand ( void ) { return -GammaPrior::rand(); }
        PsiPrior * clone ( void ) const { return new nGammaPrior(*this); }
		double mean ( void ) const { return -GammaPrior::mean(); }
};

#endif
