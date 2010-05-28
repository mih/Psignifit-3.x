/*
 *   See COPYING file distributed along with the psignifit package for
 *   the copyright and license terms
 */
#ifndef RNG_H
#define RNG_H

#include <cstdlib>
#include <cmath>
#include "errors.h"

class PsiRandom
{
	public:
		double rngcall ( void ) { return drand48(); }
		virtual double draw ( void ) { throw NotImplementedError(); }
};

class GaussRandom : public PsiRandom
{
	private:
		double mu;
		double sigma;
		bool good;
		double x1;
		double x2;
		double w;
		double y;
	public:
		GaussRandom ( double mean=0, double standarddeviation=1 ) : mu ( mean ), sigma ( standarddeviation ), good ( false ) {}
		double draw ( void );              ///< draw a random number using box muller transform
};

class UniformRandom : public PsiRandom
{
	private:
		double lower;
		double upper;
	public:
		UniformRandom ( double low=0, double up=1 ) : lower(low), upper(up) {}
		double draw ( void ) { return  (upper-lower)*rngcall()+lower; }
};

class BinomialRandom : public PsiRandom
{
	private:
		int n;
		double p;
	public:
		BinomialRandom ( int number, double probability ) : n(number), p(probability) {}
		double draw ( void );
		void setprm ( int number, double probability ) { n = number; p = probability; }
};

void set_seed(long int seedval);

#endif
