#ifndef PRIOR_H
#define PRIOR_H

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
		double pdf ( double x ) { return ( x>lower && x<upper ? height : 0 ); }
};

#endif
