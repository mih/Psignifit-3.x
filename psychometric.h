#ifndef PSYCHOMETRIC_H
#define PSYCHOMETRIC_H

#include <vector>
#include <cmath>
#include <algorithm>
#include "core.h"
#include "sigmoid.h"
#include "errors.h"
#include "prior.h"
#include "data.h"

/**
 * Standard psychometric function model
 *
 * Psi = guessingrate + (1-guessingrate-lapserate) * Sigmoid ( x )
 *
 */
class PsiPsychometric {
	private:
		int Nalternatives;
		double guessingrate;
		PsiCore * Core;
		PsiSigmoid * Sigmoid;
		std::vector<PsiPrior*> priors;
	public:
		PsiPsychometric (
			int nAFC,
			PsiCore * core,
			PsiSigmoid * sigmoid
			);    ///< Set up a psychometric function model for an nAFC task (nAFC=1 ~> yes/no)
		~PsiPsychometric ( void );
		virtual double evaluate (
			double x,
			const std::vector<double>& prm
			) const;  ///< Evaluate the psychometric function at this position
// TODO: we probably need some derivatives here
		virtual double negllikeli ( const std::vector<double>& prm, const PsiData& data ) const;   ///< negative log likelihood
		virtual double neglpost ( const std::vector<double>& prm, const PsiData& data ) const;     ///< negative log posterior
		virtual double leastfavourable ( const std::vector<double>& prm, const PsiData& data, double cut, bool threshold=true ) const; ///< derivative of log likelihood in the least favourable direction in parameter space
		virtual double deviance ( const std::vector<double>& prm, const PsiData& data ) const; ///< deviance for a given data set and parameter constellation
		const PsiCore* getCore ( void ) { return Core; }                ///< get the core of the psychometric function
		const PsiSigmoid* getSigmoid ( void ) { return Sigmoid; }       ///< get the sigmoid of the psychometric function
		void setPrior ( int index, PsiPrior* prior );                   ///< set a Prior
		int getNalternatives ( void ) const { return Nalternatives; }         ///< get the number of alternatives (1 means yes/no)
		int getNparams ( void ) const { return (Nalternatives==1 ? 4 : 3 ); } ///< get the number of free parameters of the psychometric function
		std::vector<double> getStart ( const PsiData& data ) const ;                ///< determine a starting value using logistic regression
		double getThres ( const std::vector<double>& prm, double cut ) const { return Core->inv(Sigmoid->inv(cut),prm); }  /// get the threshold at a cut between 0 and 1
		std::vector<double> getDevianceResiduals ( const std::vector<double>& prm, const PsiData& data ) const;  ///< deviance residuals for model checking
		double getRpd ( const std::vector<double>& devianceresiduals, const std::vector<double>& prm, const PsiData& data ) const;          ///< correlation between deviance residuals and predictions
		double getRkd ( const std::vector<double>& devianceresiduals, const std::vector<double>& prm ) const;          ///< correlation between deviance residuals and block sequence
};

#include "optimizer.h"

#endif
