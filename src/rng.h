/*
 *   See COPYING file distributed along with the psignifit package for
 *   the copyright and license terms
 */
#ifndef RNG_H
#define RNG_H

#include <cstdlib>
#include <cmath>
#include "errors.h"

#define NTAB 32

class PsiRandom
{
	private:
		const long IM1;
		const long IM2;
		const double AM;
		const long IMM1;
		const long IA1;
		const long IA2;
		const long IQ1;
		const long IQ2;
		const long IR1;
		const long IR2;
		const double NDIV;
		const double EPS;
		const double RNMX;
	public:
		PsiRandom ( void ) :
			IM1 ( 2147483563 ),
			IM2 ( 2147483399 ),
			AM ( 1.0/IM1 ),
			IMM1 ( IM1-1 ),
			IA1 ( 40014 ),
			IA2 ( 40692 ),
			IQ1 ( 53668 ),
			IQ2 ( 52774 ),
			IR1 ( 12211 ),
			IR2 ( 3791 ),
			NDIV ( 1 + double(IMM1)/NTAB ),
			EPS ( 1.2e-7 ),
			RNMX ( 1.0-1.2e-7 ) {}
		double rngcall ( void );
		virtual double draw ( void ) { throw NotImplementedError(); }
		virtual PsiRandom * clone ( void ) const {throw NotImplementedError(); }
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
		PsiRandom * clone ( void ) const { return new GaussRandom(*this); }
};

class UniformRandom : public PsiRandom
{
	private:
		double lower;
		double upper;
	public:
		UniformRandom ( double low=0, double up=1 ) : lower(low), upper(up) {}
		double draw ( void ) { return  (upper-lower)*rngcall()+lower; }
		PsiRandom * clone ( void ) const { return new UniformRandom(*this); }
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
		PsiRandom * clone ( void ) const { return new BinomialRandom(*this); }
};

void setSeed(long int seedval);

#endif
