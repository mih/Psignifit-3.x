/*
 *   See COPYING file distributed along with the psignifit package for
 *   the copyright and license terms
 */
#ifndef CORE_H
#define CORE_H

#include <vector>
#include "errors.h"
#include "special.h"
#include "data.h"

/** \brief inner function of the sigmoid term of the psychometric function
 *
 * The psychometric function is parameterized by two classes. The outer (PsiSigmoid) takes care
 * of the saturating nonlinearity. The PsiCore class performs some (potentially parameter dependent)
 * internal transformations of this nonlinearity.
 *
 * The PsiCore class itself is completely virtual: It is meant to be the base class off all other
 * core objects.
 */
class PsiCore
{
	public:
		virtual double g (
			double x,                       ///< stimulus intensity
			const std::vector<double>& prm  ///< parameter vector
			) { throw NotImplementedError(); }          ///< evaluate the core of the sigmoid
		virtual double dg (
			double x,                       ///< stimulus intensity
			const std::vector<double>& prm, ///< parameter vector
			int i                           ///< index of the parameter to which the derivative should be evaluated
			) { throw  NotImplementedError(); }        ///< evaluate the first derivative of the core with respect to parameter i
		virtual double ddg (
			double x,                       ///< stimulus intensity
			const std::vector<double>& prm, ///< parameter vector
			int i,                          ///< index of the first parameter to which the derivative should be evaluated
			int j                           ///< index of the second parameter to which the derivative should be evaluated
			) { throw NotImplementedError(); }         ///< evaluate the second derivative of the core with respect to parameter i and j
		virtual double inv (
			double y,                       ///< transformed intensity
			const std::vector<double>& prm  ///< parameter vector
			) { throw NotImplementedError(); }         ///< invert the core
		virtual double dinv (
			double p,                        ///< transformed inensity at which to evaluate the derivative
			const std::vector<double>& prm,  ///< parameter vector
			int i                            ///< evaluate the derivative with respect to parameter i
			) { throw NotImplementedError(); }         ///< derivative of the inverse core with respect to parameters
		virtual std::vector<double> transform (
				int nprm,                    ///< number of parameters in the final parameter vector
				double a,                    ///< intercept of the logistic regression model
				double b                     ///< slope of the logistic regression model
				) {throw NotImplementedError();}       ///< transform parameters from logistic regression to those used for this core
        virtual PsiCore * clone() const { throw NotImplementedError(); }
        static std::string getDescriptor ( void ) { throw NotImplementedError(); }///< get a short string that identifies the type of core
};

/** \brief a-b parameterization of the psychometric function
 *
 * In the original psignifit release, the nonlinearity was usually defined as a cumulative distribution function. In that
 * case two parameters describing the mean alpha and the standard deviation beta of this distribution were required. This
 * yielded a core object of the form (x-alpha)/beta. This type of internal parameterization is implemented here.
 *
 * The parameter vector is in any case expected to have the first two parameters alpha and beta
 */
class abCore : public PsiCore
{
	private:
	public:
        abCore( const PsiData* data=NULL, ///< ignored
                const int sigmoid=1,      ///< ignored
                const double alpha=0.1    ///< ignored
                ) {}                    ///< construcor
        abCore( const abCore& original) {}                    ///< copy construcor
		double g (
			double x,                        ///< stimulus intensity
			const std::vector<double>& prm   ///< parameter vector
			) { return (x-prm[0])/prm[1]; }            ///< evaluate the core of the sigmoid
		double dg (
			double x,                        ///< stimulus intensity
			const std::vector<double>& prm,  ///< parameter vector
			int i                            ///< index of the parameter to which the derivative should be evaluated
			);                                         ///< evaluate the first derivative of the core with respect to parameter i
		double ddg (
			double x,                        ///< stimulus intensity
			const std::vector<double>& prm,  ///< parameter vector
			int i,                           ///< index of the parameter to which the first derivative should be evaluated
			int j                            ///< index of the parameter to which the second derivative should be evaluated
			);                                         ///< evaluate the second derivative of the core with respect to parameters i and j
		double inv (
			 double y,                       ///< transformed intensity
			 const std::vector<double>& prm  ///< parameter vector
			 );                                        ///< invert the core
		double dinv (
			double p,                        ///< transformed intenstiy at which to evaluate the derivative
			const std::vector<double>& prm,  ///< parameter vector
			int i                            ///< evaluate the derivative with respect to parameter i
			);                                         ///< derivative of the inverse core with respect to parameter i
		std::vector<double> transform (
			int nprm,                        ///< number of parameters in the final parameter vector
			double a,                        ///< intercept of the logistic regression model
			double b                         ///< slope of the logistic regression model
			);                                         ///< transform parameters from a logistic regression model to the parameters used here
        PsiCore * clone() const;             ///< clone by value
        static std::string getDescriptor ( void ) {
            return "ab";
        }
};

/** \brief m-w parameterization of the psychmetric function
 *
 * An alternative way to parameterize the psychometric function is to describe it in terms of a threshold (m) and the width
 * of part of the function over which there is significant performance increase. What exactly "significant performance increase"
 * means is defined by a parameter alpha. By definition significant performance increase happens over the range where f(g(x|theta)) is
 * larger than alpha but smaller than 1-alpha. Obviously this definition depends on the sigmoid that is used.
 */
class mwCore : public PsiCore
{
	private:
		int sigmtype;
		double alpha;
		double zalpha;
		double zshift;
	public:
        mwCore( const PsiData* data=NULL, ///< ignored
                const int sigmoid=1,      ///< Type of the sigmoid (1=logistic, 2=gauss, 3=gumbel)
                const double alpha=0.1    ///< alpha parameter defining what "significant performance increase" means
                );                        ///< construcor
        mwCore( const mwCore& original ) :  sigmtype(original.sigmtype),
                                            alpha(original.alpha),
                                            zalpha(original.zalpha),
                                            zshift(original.zshift) {} ///< copy constructor
		double g (
			double x,                        ///< stimulus intensity
			const std::vector<double>& prm   ///< parameter vector
			);                                          ///< evaluate the core of the sigmoid
		double dg (
			double x,                        ///< stimulus intensity
			const std::vector<double>& prm,  ///< parameter vector
			int i                            ///< index of the parameter to which the derivative should be evaluated
			);                                         ///< evaluate the first derivative of the core with respect to parameter i
		double ddg (
			double x,                        ///< stimulus intensity
			const std::vector<double>& prm,  ///< parameter vector
			int i,                           ///< index of the parameter to which the first derivative should be evaluated
			int j                            ///< index of the parameter to which the second derivative should be evaluated
			);                                         ///< evaluate the second derivative of the core with respect to parameters i and j
		double inv (
			 double y,                       ///< transformed intensity
			 const std::vector<double>& prm  ///< parameter vector
			 );                                        ///< invert the core
		double dinv (
			double p,                        ///< transformed intenstiy at which to evaluate the derivative
			const std::vector<double>& prm,  ///< parameter vector
			int i                            ///< evaluate the derivative with respect to parameter i
			);                                         ///< derivative of the inverse core with respect to parameter i
		std::vector<double> transform (
			int nprm,                        ///< number of parameters in the final parameter vector
			double a,                        ///< intercept of the logistic regression model
			double b                         ///< slope of the logistic regression model
			);                                         ///< transform parameters from a logistic regression model to the parameters used here
        PsiCore * clone() const;             ///< clone by value
        static std::string getDescriptor ( void ) {
            return "mw";
        }
        double getAlpha( void ) const {
            return alpha;
        }
};

/** \brief linear core
 *
 * The core of the sigmoid is simply a*x+b, where a and b are the first two parameters. This is the parameterization that would
 * be used in the context of generalized linear models. The parameters do not have an obvious interpretation in terms of
 * psychophysically meaningful quantities. However, it might well be that in this form, the parameters are more independent, which
 * is particularly important for MCMC.
 */
class linearCore : public PsiCore
{
	public:
        linearCore( const PsiData* data=NULL, ///< ignored
                const int sigmoid=1,          ///< ignored
                const double alpha=0.1        ///< ignored
                ) {}                          ///< construcor
        linearCore( const linearCore& original ) {} ///< copy constructor
		double g (
			double x,                           ///< stimulus intensity
			const std::vector<double>& prm      ///< parameter vector
			) { return prm[0] * x + prm[1]; }   ///< evaluate the core of the sigmoid
		double dg (
			double x,                           ///< stimululs intensity
			const std::vector<double>& prm,     ///< parameter vector
			int i                               ///< index of the parameter we want the derivative to
			) { switch (i) { case 0: return x; break; case 1: return 1; break; default: return 0; break; } } ///< first derivative w.r.t. parameter i
		double ddg (
			double x,                           ///< stimulus intensity
			const std::vector<double>& prm,     ///< parameter vector
			int i,                              ///< index of the parameter we want for the first derivative
			int j                               ///< index of the parameter we want for the second derivative
			) { return 0; }                     ///< second derivative w.r.t. parameters i and j
		double inv (
			double y,                           ///< value to be inverted
			const std::vector<double>& prm      ///< parameter vector
			) { return (y-prm[1])/prm[0]; }     ///< inverse of the core
		double dinv (
			double y,                           ///< value at which the derivative of the inverse should be evaluated
			const std::vector<double>& prm,     ///< parameter vector
			int i                               ///< index of the parameter we want the derivative to
			) { switch (i) { case 0: return (prm[1]-y)/(prm[0]*prm[0]); break; case 1: return -1./prm[0]; break; default: return 0; break; } } ///< deriviative of the inverse w.r.t. parameter i
		std::vector<double> transform (
			int nprm,                           ///< number of parameters in the whole model
			double a,                           ///< intercept parameter of the logistic regression model
			double b                            ///< slope parameter of the logistic regression
			) { std::vector<double> out (nprm,0); out[0] = b; out[1] = a; return out; }   ///< transform logistic regression parameters to useful ones for this core
        PsiCore * clone() const {
            return new linearCore(*this);
        };             ///< clone by value
        static std::string getDescriptor ( void ) {
            return "linear";
        }
};

/** \brief logarithmic core
 *
 * The Weibull function typically gives a good fit for data from visual experiments. Unfortunately, the weibull distribution function
 * does not allow for a straight forward fit using generalized linear models. However, the weibull distribution function is obtained
 * if a gumbel is fit on logarithmic contrast values. This core is the same as the linearCore but for the logarithm of x
 */
class logCore : public PsiCore
{
	private:
		double scale;
	public:
        logCore( const PsiData* data=NULL, ///< use a data set to determine the correct scaling factors of initial values and initialize the objec
                const int sigmoid=1,          ///< ignored
                const double alpha=0.1        ///< ignored
                );                       ///< construcor
        logCore ( const logCore& original) : scale(original.scale) {} ///< copy constructor
		double g   (
			double x,                                 ///< stimulus intensity
			const std::vector<double>& prm            ///< parameter vector
			) throw(BadArgumentError);   ///< evaluate the core
		double dg  (
			double x,                                 ///< stimulus intensity
			const std::vector<double>& prm,           ///< parameter vector
			int i                                     ///< parameter with respect to which the derivative should evaluated
			);                   ///< evaluate derivative of the core
		double ddg (
			double x,                                 ///< stimulus intensity
			const std::vector<double>& prm,           ///< parameter vector
			int i,                                    ///< first parameter with respect to which the derivative should be taken
			int j                                     ///< second parameter with respect to which the derivative should be taken
			) { return 0; }      ///< evaluate 2nd derivative of the core
		double inv (
			double y,                                 ///< value at which to evaluate the inverse
			const std::vector<double>& prm            ///< parameter vector
			) { return exp((y-prm[1])/prm[0]); }      ///< invert the core
		double dinv (
			double y,                                 ///< value at which to evaluate the inverse
			const std::vector<double>& prm,           ///< parameter vector
			int i                                     ///< take derivative of the inverse core with respect to parameter i
			);                   ///< evaluate derivative of the inverse core with respect to parameter i
		std::vector<double> transform (
				int nprm,                             ///< number of parameters in the final model
				double a,                             ///< intercept of the logistic regression model
				double b                              ///< slope of the logistic regression model
			);                   ///< transform parameters from a logistic regression model to starting values
        PsiCore * clone() const ; ///< clone by value
        static std::string getDescriptor ( void ) {
            return "log";
        }
};

/** \brief Core for the psychofun Weibull parameterization
 *
 * The R-package psychofun by Kuss et al (2005, J Vis) uses a slightly different parameterization of the Weibull, that is parameterized in
 * terms of "threshold location m and slope at threshold s". This core should be combined with a PsiGumbelL sigmoid to obtain the Weibull or
 * with a PsiGumbelR sigmoid to obtain the reversed Weibull. However, any other combination is also valid. However, in that case, the parameters
 * m and s might not be as interpretable as they are in case of the weibull.
 */
class weibullCore : public PsiCore
{
	private:
		double twooverlog2;
		double loglog2;
		double loglina;
		double loglinb;
	public:
        weibullCore( const PsiData* data=NULL, ///< use a data set to determine the correct scaling factors of initial values and initialize the objec
                const int sigmoid=1,          ///< ignored
                const double alpha=0.1        ///< ignored
                );                       ///< construcor
		weibullCore ( const weibullCore& original ) : twooverlog2(original.twooverlog2),
                                                      loglog2(original.loglog2),
                                                      loglina(original.loglina),
                                                      loglinb(original.loglinb) {} ///< copy constructor
		double g (
			double x,                           ///< stimulus intensity
			const std::vector<double>& prm      ///< parameter vector (m,s,...)
			) { return twooverlog2*prm[0]*prm[1] * (log(x)-log(prm[0])) + loglog2; } ///< evaluate the weibull core
		double dg (
			double x,                           ///< stimulus intensity
			const std::vector<double>& prm,     ///< parameter vector
			int i                               ///< index of the parameter with respect to which the derivative should be evaluated
			) throw(BadArgumentError) ;    ///< evaluate the derivateive of the core
		double ddg (
			double x,                           ///< stimulus intenstiy
			const std::vector<double>& prm,     ///< parameter vector
			int i,                              ///< first parameter with respect to which the derivative should be taken
			int j                               ///< second parameter with respect to which the derivative should be taken
			) throw(BadArgumentError) ;            ///< evaluate the 2nd derivative of the core
		double inv (
			double y,                           ///< value at which to evaluate the inverse
			const std::vector<double>& prm      ///< parameter vector
			);           ///< invert the core
		double dinv (
			double y,                           ///< value at which to evaluate the inverse
			const std::vector<double>& prm,     ///< parameter vector
			int i                               ///< take the derivative of the inverse core with respect to parameter i
			);           ///< evaluate the derivative of the inverse core with respect to parameter i
		std::vector<double> transform (
			int nprm,                           ///< number of parameters in the final model
			double a,                           ///< intercept of the logistic regression model
			double b                            ///< slope of the logistic regression model
			);          ///< transform the parameters from a logistic regression model to starting values
        PsiCore * clone() const ; ///< clone by value
        static std::string getDescriptor ( void ) {
            return "weibull";
        }
};

/** \brief polynomial Core as used for the weibull function
 *
 * The classical weibull function is parameterized as 1-exp(-(x/alpha)^beta), this core defines the (x/alpha)^beta part in this
 * parameterization. The PsiExponential sigmoid gives the 1-exp(-.) part.
 */
class polyCore : public PsiCore
{
	private:
		double x1;
		double x2;
	public:
        polyCore( const PsiData* data=NULL, ///< use a data set to determine the correct scaling factors of initial values and initialize the objec
                const int sigmoid=1,          ///< ignored
                const double alpha=0.1        ///< ignored
                );                        ///< construcor
		polyCore ( const polyCore& original ) : x1(original.x1),
                                                x2(original.x2) {} ///< copy constructor
		double g (
			double x,                                ///< stimulus intensity
			const std::vector<double>& prm           ///< parameter vector (alpha,beta, ...)
			) { return pow( x/prm[0], prm[1] ); }    ///< evaluate the polyCore
		double dg (
			double x,                                ///< stimulus intensity
			const std::vector<double>& prm,          ///< parameter vector
			int i                                    ///< index of the parameter to which the derivative should be evaluated
			);              ///< derivative of the polyCore with respect to a parameter
		double ddg (
			double x,                                ///< stimulus intensity
			const std::vector<double>& prm,          ///< parameter vector
			int i,                                   ///< index of the first derivative parameter
			int j                                    ///< index of the 2nd derivatibe parameter
			);              ///< 2nd derivative of the polyCore object with respect to parameters
		double inv (
			double y,                                ///< value for which the core should be inverted
			const std::vector<double>& prm           ///< parameter vector
			);              ///< inverse of the core
		double dinv (
			double y,                                ///< value at which to evaluate the inverse
			const std::vector<double>& prm,          ///< parameter vector
			int i                                    ///< index of the parameter for which the derivative should be evaluated
			);              ///< derivative of the inverse core
		std::vector<double> transform (
			int nprm,                                ///< number of parameters in the final model
			double a,                                ///< intercept of the logistic regression model
			double b                                 ///< slope of the logistic regression model to starting values
			);              ///< transform the parameter from a logistic regression model to starting values
        PsiCore * clone() const ; ///< clone by value
        static std::string getDescriptor ( void ) {
            return "poly";
        }
};

#endif
