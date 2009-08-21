#ifndef DATA_H
#define DATA_H

#include <cmath>
#include <vector>
#include <iostream>
#include "errors.h"

class PsiData
{
	private:
		std::vector <double> intensities;
		std::vector <int>    Ntrials;
		std::vector <int>    Ncorrect;
		std::vector <double> Pcorrect;
		std::vector <double> logNoverK;
		int Nalternatives;
	public:
		PsiData (
			std::vector<double> x,     ///< Stimulus intensities
			std::vector<int>    N,     ///< Numbers of trials presented at the respective stimulus intensities
			std::vector<int>    k,     ///< Numbers of correct trials at the respective stimulus intensities
			int nAFC                   ///< Number of response alternatives (nAFC=1 ~> yes/no task)
			);

		PsiData (
			std::vector<double> x,     ///< Stimulus intensities
			std::vector<int>    N,     ///< Numbers of trials presented at the respective stimulus intensities
			std::vector<double> p,     ///< Fraction of correct trials at the respective stimulus intensities
			int nAFC                   ///< Number of response alternatives (nAFC=1 ~> yes/no task)
			);

		void setNcorrect ( const std::vector<int>& newNcorrect );  ///< set the number of correct responses (is probably only useful for bootstrap)
		const std::vector<double>& getIntensities ( void ) const;        ///< get the stimulus intensities
		const std::vector<int>&    getNtrials ( void ) const;            ///< get the numbers of trials at the respective stimulus intensities
		const std::vector<int>&    getNcorrect ( void ) const;           ///< get the numbers of correct trials at the respective stimulus intensities
		const std::vector<double>& getPcorrect ( void ) const;           ///< get the fraction of correct trials at the respective stimulus intensities
		double getIntensity ( int i ) const;        ///< get the stimulus intensity for block i
		int getNtrials ( int i ) const;            ///< get the numbers of trials  for block i
		int getNcorrect ( int i ) const;           ///< get the numbers of correct trials for block i
		double getPcorrect ( int i ) const;           ///< get the fraction of correct trials for block i
		int getNalternatives ( void ) const;                        ///< get the number of response alternatives (1 means yes/no task)
		int getNblocks ( void ) const { return intensities.size(); }
		double getNoverK (int i ) const;             ///< return the log of NoverK
};

#endif
