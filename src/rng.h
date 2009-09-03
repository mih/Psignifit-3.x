#ifndef RNG_H
#define RNG_H

#include <cstdlib>
#include <cmath>
#include "errors.h"

class PsiRandom
{
	public:
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
		GaussRandom ( double mean=0, double standarddeviation=1 ) : mu ( mean ), sigma ( standarddeviation ) {}
		double draw ( void );              ///< draw a random number using box muller transform
};

#endif
