#ifndef MCLIST_H
#define MCLIST_H

#include <cmath>
#include <vector>
#include <algorithm>
#include <iostream>
#include "errors.h"
#include "special.h"

/** \brief basic monte carlo samples list
 *
 * This list stores monte carlo samples and deviances, nothing else.
 */
class PsiMClist
{
	private:
		std::vector< std::vector<double> > mcestimates;
		std::vector<double> deviances;
	public:
		PsiMClist (
			int N,                      ///< number of samples to be drawn
			int nprm                    ///< number of parameters in the model that is analyzed
			) : mcestimates(nprm, std::vector<double>(N) ), deviances(N) {};   ///< Initialize the list to take N samples of nprm parameters
		PsiMClist ( const PsiMClist& mclist ) : mcestimates ( mclist.mcestimates ), deviances ( mclist.deviances ) {};   ///< copy a list of mcsamples
		~PsiMClist ( ) {} ///< destructor
		std::vector<double> getEst ( int i ) const;       ///< get a single parameter estimate at sample i
		double getEst (
			int i,                                        ///< sample index
			int prm                                       ///< parameter index
			) const;                                                           ///< get a single sample of a single parameter
		void setEst (
			int i,                                        ///< index of the sample to be set
			const std::vector<double> est,                ///< parameter vector to be set at index
			double deviance                               ///< deviance associated with the sample
			);                                                                 ///< set a sample of parameters
		virtual void setdeviance ( int i, double deviance );                   ///< set the deviance separately for sample i
		virtual double getPercentile (
			double p,                                     ///< desired percentile (in the range (0,1))
			int prm                                       ///< index of the paramter of interest
			);                                                                 ///< get a percentile for parameter prm
		virtual double getMean (
			unsigned int prm                              ///< index of the parameter of interest
			) const ;                                                          ///< get the average of parameter prm
		double getdeviance ( int i ) const;                                    ///< get the deviance of sample i
		int getNsamples ( void ) const { return mcestimates[0].size(); }       ///< get the total number of samples
		int getNparams ( void ) const { return mcestimates.size(); }           ///< get the number of parameters
		double getDeviancePercentile ( double p );                             ///< get the p-percentile of the deviance (p in the range (0,1) )
};

/** \brief list of bootstrap samples
 *
 * Bootstrap samples support some special operations that regular monte carlo samples don't. In particular bootstrap
 * samples from a psychometric function should incorporate information about
 *
 * 1. The thresholds that are associated with each parameter vector
 * 2. the bootstrap samples themselves and not only the resulting parameter estimates
 * 3. correlations of the psychometric function with the bootstrap samples "sequence"
 */
class BootstrapList : public PsiMClist
{
	private:
		bool BCa;
		std::vector<double> acceleration;
		std::vector<double> bias;
		std::vector< std::vector<int> > data;
		std::vector<double> cuts;
		std::vector< std::vector<double> > thresholds;
		std::vector<double> Rpd;
		std::vector<double> Rkd;
	public:
		BootstrapList (
			unsigned int N,                                              ///< number of samples to be drawn
			unsigned int nprm,                                           ///< number of parameters in the model
			unsigned int nblocks,                                        ///< number of blocks in the experiment
			std::vector<double> Cuts                                     ///< performance levels at which thresholds should be determined
			) : PsiMClist (N,nprm),
				BCa(false),
				acceleration(Cuts.size()),
				bias(Cuts.size()),
				data(N,std::vector<int>(nblocks)),
				cuts(Cuts),
				thresholds (Cuts.size(), std::vector<double> (N)),
				Rpd(N),
				Rkd(N)
			{ }; ///< set up the list
		// TODO: should setBCa be private and friend of parametric bootstrap?
		void setBCa (
			unsigned int i,                                               ///< index of the cut for which Bias and Acceleration should be set
			double Bias,                                                  ///< Bias to be set
			double Acceleration                                           ///< Acceleration to be set
			) { BCa=true; bias[i] = Bias; acceleration[i] = Acceleration; }  ///< set bias and acceleration to get BCa confidence intervals
		void setData (
			unsigned int i,                                               ///< index of the bootstrap sample to be set
			const std::vector<int>& newdata                                ///< responses in the new bootstrap sample
			);   ///< store a simulated data set
		std::vector<int> getData ( unsigned int i ) const;                 ///< get a simulated data set at posititon i
		double getThres ( double p, unsigned int cut );                    ///< get the p-th percentile associated with the threshold at cut
		double getThres_byPos ( unsigned int i, unsigned int cut );        ///< get the threshold for the i-th sample
		void setThres ( double thres, unsigned int i, unsigned int cut );  ///< set the value of a threshold associated with the threshold at cut
		int getNblocks ( void ) const { return data[0].size(); }           ///< get the number of blocks in the underlying dataset
		double getCut ( unsigned int i ) const;                            ///< get the value of cut i
		double getAcc ( unsigned int i ) const { return acceleration[i]; };///< get the acceleration constant for cut i
		double getBias ( unsigned int i ) const { return bias[i]; };       ///< get the bias for cut i
		// TODO: should setRpd be private and friend of parametricbootstrap?
		void setRpd ( unsigned int i, double r_pd );                       ///< set correlation between predicted values and deviance residuals for a simulated dataset
		double getRpd ( unsigned int i ) const;                            ///< get correlation between predicted values and deviance residuals for simulated dataset i
		double percRpd ( double p );                                       ///< get the p-th percentile of the correlations between predicted values and deviance residuals
		// TODO: should setRkd be private and friend of parametric bootstrap?
		void setRkd ( unsigned int i, double r_kd );                       ///< set correlation between block index and deviance residuals for a simulated dataset
		double getRkd ( unsigned int i ) const;                            ///< get correlation between block index and deviance residuals for simulated dataset i
		double percRkd ( double p );                                       ///< get the p-th percentile of the correlations between block index and deviance residuals
};

/** \brief list of JackKnife data
 *
 * JackKnifeing is not suggested for the assessment of confidence intervals or variability. Instead the close link between jackknife samples
 * and individual data points is useful to determine influential data points and outliers.
 */
class JackKnifeList : public PsiMClist
{
	private:
		double maxdeviance;
	public:
		JackKnifeList (
			unsigned int nblocks,                                             ///< number of blocks in the experiment
			unsigned int nprm,                                                ///< number of parameters in the model
			double maxlest                                                    ///< deviance of the maximum likelihood estimate on the full dataset
			) : PsiMClist ( nblocks, nprm ), maxdeviance(maxlest) {}    ///< constructor
		unsigned int getNblocks ( void ) const { return getNsamples(); } ///< get the number of blocks in the current experiment
		/** determination of influential observations is performed by checking whether a parameter changes significantly (as defined by
		 * the confidence intervals) if one observation is omitted. Thus, if leaving out one observation results in significant changes
		 * in the estimated parameters, this observation is considered "influential".
		 *
		 * \param block     index of the block to be checked
		 * \param ci_lower  lower confidence limits for each parameter in the model
		 * \param ci_upper  upper confidence limits for each parameter in the model
		 *
		 * \return true if block presents an influential observation
		 */
		bool influential ( unsigned int block, const std::vector<double>& ci_lower, const std::vector<double>& ci_upper ) const ;
		/** determination of outliers is based on the following idea: We add a new parameter that fits the data in block perfectly.
		 * If this "modified" model is significantly better than the original model, then this block is considered an outlier.
		 *
		 * \param block      index of the block to be checked
		 *
		 * \return true if block presents an outlier
		 */
		bool outlier ( unsigned int block ) const ; ///< is block an outlier?
};

#endif
