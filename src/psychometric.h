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
#include "linalg.h"

/** \brief Standard psychometric function model
 *
 * Standard model for the psychometric function that assumes that the number of correct responses is a
 * binomial random variable with parameters N (number of trials) and Psi, where Psi is
 *
 * Psi = guessingrate + (1-guessingrate-lapserate) * Sigmoid ( x | theta )
 *
 * For an nAFC task, the guessingrate is typicall fixed at 1/n.
 *
 * The term Sigmoid ( x | theta ) is represented by two objects:
 * a PsiSigmoid f, that describes a nonlinear function from the real numbers to (0,1) and
 * a PsiCore g, that describes the "internal" workings of the nonlinear function. Thus,
 * the Term Sigmoid ( x | theta ) equals f(g(x,theta)).
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
			int nAFC,                                                                ///< number of alternatives in the task (1 indicating yes/no)
			PsiCore * core,                                                          ///< internal part of the nonlinear function (in many cases this is actually a linear function)
			PsiSigmoid * sigmoid                                                     ///< "external" saturating part of the nonlinear function
			);    ///< Set up a psychometric function model for an nAFC task (nAFC=1 ~> yes/no)
		~PsiPsychometric ( void ); ///< destructor
		virtual double evaluate (
			double x,                                                                ///< stimulus intensity
			const std::vector<double>& prm                                           ///< parameters of the psychometric function model
			) const;  ///< Evaluate the psychometric function at this position
		virtual double negllikeli (
			const std::vector<double>& prm,                                          ///< parameters of the psychometric function model
			const PsiData* data                                                      ///< data for which the likelihood should be evaluated
			) const;   ///< negative log likelihood
		virtual double neglpost (
			const std::vector<double>& prm,                                          ///< parameters of the psychometric function model
			const PsiData* data                                                      ///< data for which the posterior should be evaluated
			) const;     ///< negative log posterior  (unnormalized)
		virtual double leastfavourable (
			const std::vector<double>& prm,                                          ///< parameters of the psychometric function model
			const PsiData* data,                                                     ///< data for which the likelihood should be evaluated
			double cut,                                                              ///< performance level at which the threshold should be evaluated
			bool threshold=true                                                      ///< should the calculations be performed for thresholds? (anything else is not yet implemented)
			) const; ///< derivative of log likelihood in the least favourable direction in parameter space
		virtual double deviance (
			const std::vector<double>& prm,                                          ///< parameters of the psychometric functin model
			const PsiData* data                                                      ///< data for which the likelihood should be evaluated
			) const; ///< deviance for a given data set and parameter constellation
		virtual Matrix * ddnegllikeli (
				const std::vector<double>& prm,                                      ///< parameters at which the second derivative should be evaluated
				const PsiData* data                                                  ///< data for which the likelihood should be evaluated
				) const;                                          ///< 2nd derivative of the negative log likelihood (newly allocated matrix)
		virtual std::vector<double> dnegllikeli (
				const std::vector<double>& prm,                                      ///< parameters at which the first derivative should be evaluated
				const PsiData* data                                                  ///< data for which the likelihood should be evaluated
				) const;                                          ///< 1st derivative of the negative log likelihood
		const PsiCore* getCore ( void ) { return Core; }                ///< get the core of the psychometric function
		const PsiSigmoid* getSigmoid ( void ) { return Sigmoid; }       ///< get the sigmoid of the psychometric function
		void setPrior ( int index, PsiPrior* prior );                   ///< set a Prior for the parameter indicated by index
		int getNalternatives ( void ) const { return Nalternatives; }         ///< get the number of alternatives (1 means yes/no)
		int getNparams ( void ) const { return (Nalternatives==1 ? 4 : 3 ); } ///< get the number of free parameters of the psychometric function
		std::vector<double> getStart ( const PsiData* data ) const ;                ///< determine a starting value using logistic regression on a dataset
		double getThres (
			const std::vector<double>& prm,                                          ///< parameters of the psychometric function model
			double cut                                                               ///< performance level at which the threshold should be evaluated
			) const { return Core->inv(Sigmoid->inv(cut),prm); }  ///< get the threshold at a cut between 0 and 1
		std::vector<double> getDevianceResiduals (
			const std::vector<double>& prm,                                          ///< parameters of the psychometric function model
			const PsiData* data                                                      ///< data for which the deviance residuals should be determined
			) const;  ///< deviance residuals for model checking
		double getRpd (
			const std::vector<double>& devianceresiduals,                            ///< deviance residuals as determined by getDevianceResiduals
			const std::vector<double>& prm,                                          ///< parameters of the psychometric function model
			const PsiData* data                                                      ///< data set corresponding to the deviance residuals
			) const;          ///< correlation between deviance residuals and predictions
		double getRkd ( const std::vector<double>& devianceresiduals ) const;        ///< correlation between deviance residuals and block sequence
		double dllikeli (
			std::vector<double> prm,                                                     ///< parameters of the model
			const PsiData* data,                                                         ///< data for which the likelihood should be evaluated
			unsigned int i                                                               ///< index of the parameter for which the derivative should be evaluated
			) const;                                                                 ///< derivative of the negative loglikelihood with respect to parameter i
		double dlposteri (
			std::vector<double> prm,                                                     ///< parameters of the psychometric function model
			const PsiData* data,                                                         ///< data for which the likelihood should be valuated
			unsigned int i                                                               ///< index of the parameter for which the derivative should be evaluated
			) const;                                                                 ///< derivative of the negative log posterior with respect to parameter i
};

#include "optimizer.h"

#endif
