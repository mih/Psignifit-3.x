#ifndef CORE_H
#define CORE_H

#include <vector>
#include "errors.h"
#include "special.h"

class PsiCore
{
	public:
		virtual double g (
			double x,                 ///< stimulus intensity
			const std::vector<double>& prm ///< parameter vector
			) { throw NotImplementedError(); }          ///< evaluate the core of the sigmoid
		virtual double dg (
			double x,                  ///< stimulus intensity
			const std::vector<double>& prm, ///< parameter vector
			int i                      ///< index of the parameter to which the derivative should be evaluated
			) { throw  NotImplementedError(); } ///< evaluate the first derivative of the core with respect to parameter i
		virtual double ddg (
			double x,                  ///< stimulus intensity
			const std::vector<double>& prm, ///< parameter vector
			int i,                     ///< index of the first parameter to which the derivative should be evaluated
			int j                      ///< index of the second parameter to which the derivative should be evaluated
			) { throw NotImplementedError(); } ///<evaluate the second derivative of the core with respect to parameter i and j
		virtual double inv (
			double y,                  ///< transformed intensity
			const std::vector<double>& prm  ///< parameter vector
			) { throw NotImplementedError(); } ///< invert the core
		virtual double dinv (
			double p,                   ///< transformed inensity at which to evaluate the derivative
			const std::vector<double>& prm,  ///< parameter vector
			int i                       ///< evaluate the derivative with respect to parameter i
			) { throw NotImplementedError(); } ///< derivative of the inverse core with respect to parameters
		virtual std::vector<double> transform ( int nprm, double a, double b ) {throw NotImplementedError();} ///< transform parameters from logistic regression to those used for this core
};

class abCore : public PsiCore
{
	private:
	public:
		double g    ( double x, const std::vector<double>& prm ) { return (x-prm[0])/prm[1]; }
		double dg   ( double x, const std::vector<double>& prm, int i);
		double ddg  ( double x, const std::vector<double>& prm, int i, int j);
		double inv  ( double y, const std::vector<double>& prm );
		double dinv ( double p, const std::vector<double>& prm, int i );
		std::vector<double> transform ( int nprm, double a, double b );
};

class mwCore : public PsiCore
{
	private:
		int sigmtype;
		double alpha;
		double zalpha;
		double zshift;
	public:
		mwCore ( int sigmoid, double al=0.1 );
		double g ( double x, const std::vector<double>& prm );
		double dg ( double x, const std::vector<double>& prm, int i );
		double ddg ( double x, const std::vector<double>& prm, int i, int j );
		double inv ( double y, const std::vector<double>& prm );
		double dinv ( double p, const std::vector<double>& prm, int i );
		std::vector<double> transform ( int nprm, double a, double b );
};

#endif
