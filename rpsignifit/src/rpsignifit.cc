#include "psipp.h"
#include <R.h>
#include <Rmath.h>
#include <vector>
#include <cstdio>

PsiSigmoid * determine_sigmoid ( char *sigmoid ) {
	if ( !strcmp(sigmoid,"logistic") ) {
		return new PsiLogistic ();
	} else if ( !strcmp(sigmoid,"exp") ) {
		return new PsiExponential ();
	} else if ( !strcmp(sigmoid,"gauss") ) {
		return new PsiGauss ();
	} else if ( !strcmp(sigmoid,"gumbel_l") || !strcmp(sigmoid,"gumbel") || !strcmp(sigmoid,"lgumbel") ) {
		return new PsiGumbelL ();
	} else if ( !strcmp(sigmoid,"gumbel_r") || !strcmp(sigmoid,"rgumbel") ) {
		return new PsiGumbelL();
	} else if ( !strcmp(sigmoid,"cauchy") ) {
		return new PsiCauchy();
	} else {
		Rprintf ( "WARNING: no valid sigmoid!" );
		return NULL;
	}
}

PsiCore * determine_core ( char *core, const PsiSigmoid* Sigmoid, const PsiData* data ) {
	double dummy;
	if ( !strcmp(core,"ab") ) {
		return  new abCore();
	} else if ( !strncmp(core,"mw",2) ) {
		sscanf(core, "mw%lf", &dummy);
		return  new mwCore(Sigmoid->getcode(), dummy);
	} else if ( !strcmp(core,"linear") ) {
		return  new linearCore();
	} else if ( !strcmp(core,"log") ) {
		return  new logCore(data);
	} else if ( !strcmp(core,"poly") ) {
		return  new polyCore(data);
	} else if ( !strcmp(core,"weibull") ) {
		return  new weibullCore(data);
	} else {
		Rprintf ( "WARNING: no valid core!" );
		return NULL;
	}
}

PsiPrior * determine_prior ( char *prior ) {
	double prm1,prm2;
	if ( !strncmp(prior,"Beta",4) ) {
		sscanf(prior,"Beta(%lf,%lf)",&prm1,&prm2);
		return new BetaPrior(prm1,prm2);
	} else if ( !strncmp(prior,"Gamma",5) ) {
		sscanf(prior,"Gamma(%lf,%lf)",&prm1,&prm2);
		return new GammaPrior(prm1,prm2);
	} else if ( !strncmp(prior,"Gauss",5) ) {
		sscanf(prior,"Gauss(%lf,%lf)",&prm1,&prm2);
		return new GaussPrior(prm1,prm2);
	} else if ( !strncmp(prior,"Uniform",7) ) {
		sscanf(prior,"Uniform(%lf,%lf)",&prm1,&prm2);
		return new UniformPrior(prm1,prm2);
	} else {
		return new PsiPrior();
	}
}

extern "C" {
////////////////////////////////////
// Functions go here
////////////////////////////////////

void mapestimate (
		double * x,    // stimulus intensities
		int *k,        // response counts of correct- (nAFC) or Yes-responses (Yes/No)
		int *n,        // numbers of trials per block
		int *K,        // number of blocks
		char **sigmoid,// the sigmoid to be used
		char **core,   // core description
		int *nafc,     // number of alternatives in the task (a value < 2 indicates Yes/No)
		double *estimate, // output array for the estimated values
		int *nparams,     // number of parameters
		double *deviance, // output: deviance
		char **priors
		) {
	std::vector<double> stimulus_intensities ( *K );
	std::vector<int>    number_of_trials     ( *K );
	std::vector<int>    number_of_correct    ( *K );
	int i;

	for ( i=0; i<*K; i++ ) {
		stimulus_intensities[i] = x[i];
		number_of_trials[i]     = n[i];
		number_of_correct[i]    = k[i];
	}

	PsiData *data = new PsiData ( stimulus_intensities, number_of_trials, number_of_correct, *nafc );

	PsiSigmoid *Sigmoid = determine_sigmoid ( *sigmoid );
	if (Sigmoid==NULL) { delete data; return; }

	PsiCore *Core = determine_core ( * core, Sigmoid, data );
	if (Core==NULL) { delete data; delete Sigmoid; return; }

	PsiPsychometric * pmf = new PsiPsychometric ( *nafc, Core, Sigmoid );
	if (*nparams != pmf->getNparams() ) {
		Rprintf ( "WARNING: output vector length does not match number of parameters!" );
		delete data;
		delete Sigmoid;
		delete Core;
		delete pmf;
		return;
	}

	for ( i=0; i<*nparams; i++ ) {
		pmf->setPrior ( i, determine_prior ( priors[i] ) );
	}

	PsiOptimizer *opt = new PsiOptimizer ( pmf, data );
	std::vector<double> est = opt->optimize ( pmf, data );
	delete opt;


	for ( i=0; i<*nparams; i++ )
		estimate[i] = est[i];

	*deviance = pmf->deviance(est,data);

	delete data;
	delete pmf;

	return;
}

/*
void bootstrap ( ... ) {
}

void mcmc ( ... ) {
}
*/

// Some more?
}
