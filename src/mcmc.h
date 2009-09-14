#ifndef MCMC_H
#define MCMC_H

#include <vector>
#include "psychometric.h"
#include "rng.h"
#include "mclist.h"

class PsiSampler
{
	private:
		const PsiPsychometric * model;
		const PsiData * data;
	public:
		PsiSampler ( const PsiPsychometric * Model, const PsiData * Data ) : model(Model), data(Data) {}///< set up a sampler to sample from the posterior of the parameters of pmf given the data dat
		virtual std::vector<double> draw ( void ) { throw NotImplementedError(); }                     ///< draw a sample from the posterior
		virtual void setTheta ( const std::vector<double> theta ) { throw NotImplementedError(); }     ///< set the "state" of the underlying markov chain
		virtual std::vector<double> getTheta ( void ) { throw NotImplementedError(); }                 ///< get the "state" of the underlying markov chain
		virtual void setstepsize ( double size, unsigned int param ) { throw NotImplementedError(); }  ///< set the size of the steps for parameter param of the sampler
		virtual void setstepsize ( const std::vector<double>& sizes ) { throw NotImplementedError(); } ///< set all stepsizes of the sampler
		virtual double getDeviance ( void ) { throw NotImplementedError(); }                           ///< return the model deviance for the current state
		virtual PsiMClist sample ( unsigned int N ) { throw NotImplementedError(); }                   ///< draw N samples from the posterior
		const PsiPsychometric * getModel() const { return model; }                                     ///< return the underlying model instance
		const PsiData         * getData()  const { return data;  }                                     ///< return the underlying data instance
};

class MetropolisHastings : public PsiSampler
{
	private:
		PsiRandom* propose;
		std::vector<double> currenttheta;
		std::vector<double> newtheta;
		std::vector<double> stepwidths;
		double qold;
		double currentdeviance;
		int accept;
	public:
		MetropolisHastings (
			const PsiPsychometric * Model,                                                  ///< psychometric funciton model to sample from
			const PsiData * Data,                                                           ///< data to base inference on
			PsiRandom* proposal                                                             ///< proposal distribution (will usually be a gaussian)
			);                                                          ///< initialize the sampler
		~MetropolisHastings ( void ) { delete propose; }
		std::vector<double> draw ( void );                                                ///< perform a metropolis hastings step and draw a sample from the posterior
		void setTheta ( const std::vector<double>& prm );                                 ///< set the current state of the sampler
		std::vector<double> getTheta ( void ) { return currenttheta; }                    ///< get the current state of the sampler
		double getDeviance ( void ) { return currentdeviance; }                           ///< get the current deviance
		void setstepsize ( double size, unsigned int param );                             ///< set the standard deviation of the proposal distribution for parameter param
		void setstepsize ( const std::vector<double>& sizes );                            ///< set standard deviations of the proposal distribution for all parameters at once
		PsiMClist sample ( unsigned int N );                                              ///< draw N samples from the posterior
		unsigned int getNparams ( void ) { return newtheta.size(); }                      ///< get the number of parameters for which the sampler is set up
};

class HybridMCMC : public PsiSampler
{
	private:
		PsiRandom* proposal;
		std::vector<double> currenttheta;
		std::vector<double> newtheta;
		std::vector<double> momentum;
		double currentH;
		double newH;
		double energy;
		double newenergy;
		std::vector<double> gradient;
		std::vector<double> currentgradient;
		std::vector<double> stepsizes;
		int Nleapfrog;
		int Naccepted;
		void leapfrog ( void );
	public:
		HybridMCMC (
			const PsiPsychometric * Model,                                                  ///< psychometric function model to sample from
			const PsiData * Data,                                                          ///< data to base inference on
			int Nleap                                                                     ///< number of leapfrog steps to be performed for each sample
			);                                                             ///< initialize the sampler
		~HybridMCMC ( void ) { delete proposal; }
		std::vector<double> draw ( void );                                                ///< draw a sample from the posterior
		void setTheta ( const std::vector<double>& prm );                                 ///< set the current state of the sampler
		std::vector<double> getTheta ( void ) { return currenttheta; }                    ///< get the current state of the sampler
		void setstepsize ( double size, unsigned int param );                             ///< set stepsize of the leapfrog integration for parameter param
		void setstepsize ( const std::vector<double>& sizes );                            ///< set all stepsizes of leapfrog integration for all parameters at once
		double getDeviance ( void );                                                      ///< get the current deviance
		PsiMClist sample ( unsigned int N );                                              ///< draw N samples from the posterior
};

#endif
