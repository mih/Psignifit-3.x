#ifndef MCMC_H
#define MCMC_H

#include <vector>
#include "psychometric.h"
#include "rng.h"

class PsiMetropolisHastings
{
	private:
		const PsiPsychometric model;
		const PsiData data;
		const PsiRandom* propose;
		std::vector<double> currenttheta;
		std::vector<double> newtheta;
		std::vector<double> stepwidths;
		double qold;
		double currentdeviance;
	public:
		// TODO: Should Sampler arguments be const? Should we store copies?
		PsiMetropolisHastings ( const PsiPsychometric * pmf, const PsiData * dat );       ///< initialize the sampler
		std::vector<double> sample ( void );                                              ///< perform a metropolis hastings step and draw a sample from the posterior
		void setPrm ( std::vector<double>& prm );                                         ///< set the current state of the sampler
		std::vector<double> getPrm ( void );                                              ///< get the current state of the sampler
		double getDeviance ( void );                                                      ///< get the current deviance
};

#endif
