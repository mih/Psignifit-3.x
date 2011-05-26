#ifndef INTEGRATE_H
#define INTEGRATE_H

#include "psychometric.h"
#include "data.h"
#include "getstart.h"
#include "bootstrap.h"
#include <vector>
#include <algorithm>

// #define DEBUG_INTEGRATE

std::vector<double> raw_grid (
		const PsiData *data,
		const PsiPsychometric *pmf,
		unsigned int prmindex,
		unsigned int gridsize=7
		);

std::vector<double> cdf_grid (
		const PsiPrior *fitted_dist, ///< fitted distribution
		double pmin=0.1,       ///< minimum cumulative probability
		double pmax=0.9,       ///< maximum cumulative probability
		unsigned int gridsize=7///< number of gridpoints
		);

std::vector<double> lingrid ( double xmin, double xmax, unsigned int gridsize );

void normalize_margin ( std::vector<double>* margin );

double sd ( const std::vector<double>& x );

void integrate_3d (
		const PsiPsychometric *pmf,        ///< psychometric function model to be evaluated
		const PsiData *data,               ///< data set to be evaluated
		const std::vector<double>& grid1,  ///< grid points for first parameter
		const std::vector<double>& grid2,  ///< grid points for second parameter
		const std::vector<double>& grid3,  ///< grid points for third parameter
		std::vector<double> *margin1,      ///< evaluated marginal for first parameter
		std::vector<double> *margin2,      ///< evaluated marginal for second parameter
		std::vector<double> *margin3       ///< evaluated marginal for third parameter
		);  ///< integration on a three dimensional grid

void integrate_4d (
		const PsiPsychometric *pmf,        ///< psychometric function model to be evaluated
		const PsiData *data,               ///< data set to be evaluated
		const std::vector<double>& grid1,  ///< grid points for first parameter
		const std::vector<double>& grid2,  ///< grid points for second parameter
		const std::vector<double>& grid3,  ///< grid points for third parameter
		const std::vector<double>& grid4,  ///< grid points for fourth parameter
		std::vector<double> *margin1,      ///< evaluated marginal for first parameter
		std::vector<double> *margin2,      ///< evaluated marginal for second parameter
		std::vector<double> *margin3,      ///< evaluated marginal for third parameter
		std::vector<double> *margin4       ///< evaluated marginal for fourth parameter
		);  ///< integration on a four dimensional grid

void integrate_grid (
		const PsiPsychometric *pmf,        ///< psychometric function model to be evaluated
		const PsiData *data,               ///< data set to be evaluated
		const std::vector< std::vector<double> >& grid,   ///< all grids to be evaluated
		std::vector<double> *margin1,      ///< evaluated marginal for first parameter
		std::vector<double> *margin2,      ///< evaluated marginal for second parameter
		std::vector<double> *margin3,      ///< evaluated marginal for third parameter
		std::vector<double> *margin4=NULL  ///< evaluated marginal for fourth parameter
		) throw (PsiError);  ///< integration on a three or four dimensional grid (selects the appropriate integration function)

std::vector<double> fit_posterior (
		const std::vector<double>& x,
		const std::vector<double>& fx,
		const std::vector<double>& start,
		unsigned int index
		);

double error_gauss (
		const std::vector<double>& prm,
		const std::vector<double>& x,
		const std::vector<double>& fx
		);

double error_gamma (
		const std::vector<double>& prm,
		const std::vector<double>& x,
		const std::vector<double>& fx
		);

double error_beta (
		const std::vector<double>& prm,
		const std::vector<double>& x,
		const std::vector<double>& fx
		);

class PsiIndependentPosterior {
	private:
		unsigned int nparams;
		std::vector<PsiPrior*> fitted_posteriors;
		std::vector< std::vector<double> > grids;
		std::vector< std::vector<double> > margins;
	public:
		PsiIndependentPosterior (
				unsigned int nprm,
				std::vector<PsiPrior*> posteriors,
				std::vector< std::vector<double> > x,
				std::vector< std::vector<double> > fx
				);
			// { unsigned int i; for (i=0; i<nparams; i++) fitted_posteriors[i] = new PsiPrior; }
		PsiIndependentPosterior ( const PsiIndependentPosterior& post ) :
			nparams ( post.nparams ), fitted_posteriors ( post.nparams ), grids ( post.grids ), margins ( post.margins )
			{ unsigned int i; for ( i=0; i<nparams; i++ ) fitted_posteriors[i] = post.fitted_posteriors[i]->clone(); }
		~PsiIndependentPosterior ( void ) { unsigned int i; for ( i=0; i<nparams; i++ ) delete fitted_posteriors[i]; }
		PsiPrior *get_posterior ( unsigned int parameter ) { return fitted_posteriors[parameter]->clone(); }
		std::vector<double> get_grid ( unsigned int parameter ) { return grids[parameter]; }
		std::vector<double> get_margin ( unsigned int parameter ) { return margins[parameter]; }
};

PsiIndependentPosterior independent_marginals (
		const PsiPsychometric *pmf,
		const PsiData *data,
		unsigned int nrefinements=3,
		unsigned int gridsize=7
		);

MCMCList sample_posterior (
		const PsiPsychometric *pmf,
		const PsiData *data,
		PsiIndependentPosterior& post,
		unsigned int nsamples=600,
		unsigned int propose=25
		);

void sample_diagnostics (
		const PsiPsychometric *pmf,
		const PsiData *data,
		MCMCList *samples
		);

#endif
