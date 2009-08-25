#ifndef DATA_H
#define DATA_H

#include <cmath>
#include <vector>
#include <iostream>
#include "errors.h"

/** \brief basic data set class
 *
 * This class holds data from a psychophysical experiment. The experiment can in principle be performed as an nAFC task
 * or as a yes/no task. Typically for an nAFC task, the number of correct responses would be recorded and for a yes/no
 * task the number of yes responses would be recorded. The nomenclature of the class corresponds to the nAFC task.
 *
 * The block order and number can be crucial under some circumstances, so don't lump multiple blocks into one just
 * because they have the same stimulus intensity.
 */
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
			);                                   ///< constructor
		PsiData (
			std::vector<double> x,     ///< Stimulus intensities
			std::vector<int>    N,     ///< Numbers of trials presented at the respective stimulus intensities
			std::vector<double> p,     ///< Fraction of correct trials at the respective stimulus intensities
			int nAFC                   ///< Number of response alternatives (nAFC=1 ~> yes/no task)
			);                                   ///< constructor

		void setNcorrect (
			const std::vector<int>& newNcorrect  ///< new number of correct responses
			);  ///< set the number of correct responses (is probably only useful for bootstrap)
		const std::vector<double>& getIntensities ( void ) const;        ///< get the stimulus intensities
		const std::vector<int>&    getNtrials ( void ) const;            ///< get the numbers of trials at the respective stimulus intensities
		const std::vector<int>&    getNcorrect ( void ) const;           ///< get the numbers of correct trials at the respective stimulus intensities
		const std::vector<double>& getPcorrect ( void ) const;           ///< get the fraction of correct trials at the respective stimulus intensities
		double getIntensity ( int i ) const;                             ///< get the stimulus intensity for block i
		int getNtrials ( int i ) const;                                  ///< get the numbers of trials  for block i
		int getNcorrect ( int i ) const;                                 ///< get the numbers of correct trials for block i
		double getPcorrect ( int i ) const;                              ///< get the fraction of correct trials for block i
		int getNalternatives ( void ) const;                             ///< get the number of response alternatives (1 means yes/no task)
		int getNblocks ( void ) const { return intensities.size(); }     ///< get the number of blocks in the data set
		double getNoverK (int i ) const;                                 ///< return the log of NoverK for block i
};

#endif
