#ifndef MCLIST_H
#define MCLIST_H

#include <cmath>
#include <vector>
#include <algorithm>
#include "errors.h"
#include <iostream>
#include "special.h"

class PsiMClist
{
	private:
		std::vector< std::vector<double> > mcestimates;
		std::vector<double> deviances;
	public:
		PsiMClist ( int N, int nprm )
			: mcestimates(nprm, std::vector<double>(N) ),
				deviances(N)
			{ };   ///< Initialize the list to take N samples of nprm parameters
		PsiMClist ( ) {}
		std::vector<double> getEst ( int i ) const;                                                 ///< get a single parameter estimate
		double getEst ( int i, int prm ) const;                                                           ///< get a single sample of a single parameter
		void setEst ( int i, const std::vector<double> est, double deviance );                                             ///< set a sample of parameters
		virtual void setdeviance ( int i, double deviance );                         ///< set the deviance separately
		virtual double getPercentile ( double p, int prm );                                               ///< get a percentile for parameter prm
		double getdeviance ( int i ) const; ///< get the log likelihood of sample i
		int getNsamples ( void ) const { return mcestimates[0].size(); }
		int getNparams ( void ) const { return mcestimates.size(); }
		double getDeviancePercentile ( double p );
};

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
		BootstrapList ( unsigned int N, unsigned int nprm, unsigned int nblocks, std::vector<double> Cuts )
			: PsiMClist (N,nprm),
				BCa(false),
				acceleration(Cuts.size()),
				bias(Cuts.size()),
				data(N,std::vector<int>(nblocks)),
				cuts(Cuts),
				thresholds (Cuts.size(), std::vector<double> (N)),
				Rpd(N),
				Rkd(N)
			{ }; ///< set up the list
		void setBCa ( unsigned int i, double Bias, double Acceleration ) { BCa=true; bias[i] = Bias; acceleration[i] = Acceleration; }  ///< set bias and acceleration to get BCa confidence intervals
		void setData ( unsigned int i, const std::vector<int> newdata );   ///< store a simulated data set
		std::vector<int> getData ( unsigned int i ) const;                 ///< get a simulated data set
		double getThres ( double p, unsigned int cut );                    ///< get the p-th percentile associated with the threshold at cut
		void setThres ( double thres, unsigned int i, unsigned int cut );  ///< set the value of a threshold associated with the threshold at cut
		int getNblocks ( void ) const { return data[0].size(); }           ///< get the number of blocks in the underlying dataset
		double getCut ( unsigned int i ) const;                            ///< get the value of a cut
		double getAcc ( unsigned int i ) const { return acceleration[i]; };///< get the acceleration constant
		double getBias ( unsigned int i ) const { return bias[i]; };       ///< get the bias
		void setRpd ( unsigned int i, double r_pd );                       ///< set correlation between predicted values and deviance residuals for a simulated dataset
		double getRpd ( unsigned int i ) const;                            ///< get correlation between predicted values and deviance residuals for simulated dataset i
		double percRpd ( double p );                                       ///< get the p-th percentile of the correlations between predicted values and deviance residuals
		void setRkd ( unsigned int i, double r_kd );                       ///< set correlation between block index and deviance residuals for a simulated dataset
		double getRkd ( unsigned int i ) const;                            ///< get correlation between block index and deviance residuals for simulated dataset i
		double percRkd ( double p );                                       ///< get the p-th percentile of the correlations between block index and deviance residuals
};

class JackKnifeList : public PsiMClist
{
	private:
		double maxdeviance;
	public:
		JackKnifeList ( unsigned int nblocks, unsigned int nprm, double maxlest ) : PsiMClist ( nblocks, nprm ), maxdeviance(maxlest) {}
		unsigned int getNblocks ( void ) const { return getNsamples(); }
		bool influential ( unsigned int block, const std::vector<double>& ci_lower, const std::vector<double>& ci_upper ) const ;            ///< is block an influential observation?
		bool outlier ( unsigned int block ) const ; ///< is block an outlier?
};

#endif
