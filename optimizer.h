#ifndef OPTIMIZER_H
#define OPTIMIZER_H

#include <vector>
#include "psychometric.h"
#include "data.h"

class PsiOptimizer
{
	private:
		// some variables for the internal processing of the optimization process
		int nparameters;                             // Number of parameters
		std::vector< std::vector<double> > simplex;  // data of the simplex
		std::vector<double> fx;                      // function values at the simplex nodes
		std::vector<double> x;                       // a single simplex node
		std::vector<double> xx;                      // another single simplex node
		std::vector<double> start;                   // starting values
		std::vector<bool>   modified;                // bookkeeping vector to indicate which simplex nodes have changed, i.e. which function values need to be updated
	public:
		PsiOptimizer ( const PsiPsychometric * model, const PsiData * data ); ///< set up everything
		~PsiOptimizer ( void );                                   ///< clean up everything
		std::vector<double> optimize ( const PsiPsychometric * model, const PsiData * data ); ///< Start the optimization process
};

#endif
