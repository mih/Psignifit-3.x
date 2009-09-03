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
		mwCore (
			int sigmoid,                     ///< Type of the sigmoid (1=logistic, 2=gauss, 3=gumbel)
			double al=0.1                    ///< alpha parameter defining what "significant performance increase" means
			);                                          ///< constructor
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
};

/** \brief logarithmic core
 *
 * The Weibull function typically gives a good fit for data from visual experiments. Unfortunately, the weibull distribution function
 * does not allow for a straight forward fit using generalized linear models. However, the weibull distribution function is obtained
 * if a gumbel is fit on logarithmic contrast values. This core is the same as the linearCore but for the logarithm of x
 */
class logarithmicCore : public PsiCore
{
	private:
		double scale;
	public:
		logarithmicCore ( const PsiData* data );      ///< use a data set to determine the correct scaling factors of initial values and initialize the object
		double g   (
			double x,                                 ///< stimulus intensity
			const std::vector<double>& prm            ///< parameter vector
			) { return prm[0] * log(x) + prm[1]; }    ///< evaluate the core
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
};

#endif
