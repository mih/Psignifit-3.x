#include "bootstrap.h"
#include "rng.h"

#ifdef DEBUG_BOOTSTRAP
#include <iostream>
#endif

void newsample ( const PsiData * data, const std::vector<double>& p, std::vector<int> * sample ) {
	/* Draw a new sample from the psychometric function */
	BinomialRandom binomial ( 10, 0.5 );    // Initialize with nonsense parameters
	int k;                                            // Block index

	for ( k=0; k<data->getNblocks(); k++ ) {
		binomial.setprm ( data->getNtrials(k), p[k] );
		(*sample)[k] = binomial.draw ();
	}
}

void determineBCa ( const std::vector<double>& l_LF, const std::vector<double>& u_i, double initialthreshold, double *bias, double*acc ) {
	// Calculate BCa constants
	double E_l(0);     // expected likelihood
	double E_l3(0);    // expectation of cubic likelihood
	double var_l(0);   // variance / standard deviation of likelihood
	double w(0);       // probability to obtain thresholds less than ML-fit

	int b,B ( l_LF.size() );

	// Get expectations
	for ( b=0; b<B; b++ ) {
		E_l += l_LF[b];
		E_l3 += l_LF[b]*l_LF[b]*l_LF[b];
		w += int( u_i[b] < initialthreshold );
	}

	// Get variance/standard deviation
	for ( b=0; b<B; b++ ) {
		var_l += (l_LF[b]-E_l)*(l_LF[b]-E_l);
	}
	var_l /= B-1;
	var_l = sqrt(var_l); // note that this is not variance anymore but standard deviation

#ifdef DEBUG_BOOTSTRAP
	std::cerr << " E_l3=" << E_l3 << "  var_l="<< var_l << "\n";
	std::cerr << " w = " << w/(B+1) << " real_w = " << invPhi(w/(B+1)) << "\n";
#endif

	// output
	*bias = invPhi(w/(B+1));
	*acc  = E_l3 / (6*var_l*var_l*var_l);
}

BootstrapList parametricbootstrap ( int B, const PsiData * data, const PsiPsychometric* model, std::vector<double> cuts, std::vector<double>* param, bool BCa )
{
#ifdef DEBUG_BOOTSTRAP
	std::cerr << "Starting bootstrap\n Cuts size=" << cuts.size() << " "; std::cerr.flush();
#endif
	BootstrapList bootstrapsamples ( B, model->getNparams(), data->getNblocks(), cuts );
	int b,k,l,cut;                               // iteration variables for bootstrap sample, block, l-general purpose third level iteration, cut
	std::vector< std::vector<double> > l_LF (cuts.size(), std::vector<double>(B));   // vector of double-vectors
	std::vector< std::vector<double> > u_i  (cuts.size(), std::vector<double>(B));
	PsiOptimizer opt ( model, data );                          // for ML-Fitting
	PsiData * localdataset = new PsiData ( data->getIntensities(),  // local because it changes in every iteration
			data->getNtrials(),
			data->getNcorrect(),
			data->getNalternatives() );

	std::vector<double> initialfit ( model->getNparams() );       // generating parameters for the bootstrap samples
	if (param==NULL) {
		initialfit = opt.optimize( model, data );
	} else
		initialfit = *param;
	std::vector<double> p          ( data->getNblocks() );       // predicted p-correct for parametric bootstrap
	for ( k=0; k<data->getNblocks(); k++ ) { p[k] = model->evaluate ( data->getIntensity(k), initialfit ); }

	std::vector<double> localfit   ( model->getNparams() );
	std::vector<int>    sample     ( data->getNblocks() );
	std::vector<double> initialthresholds ( cuts.size() );
	std::vector<double> devianceresiduals ( data->getNblocks() );
	double deviance;

	for (cut=0; cut<cuts.size(); cut++)
		initialthresholds[cut] = model->getThres(initialfit,cuts[cut]);

	for ( b=0; b<B; b++ ) {
		// Resampling
		newsample ( data, p, &sample );         // draw a new sample
		localdataset->setNcorrect ( sample );   // put the new sample to the localdataset
		bootstrapsamples.setData ( b, sample ); // store the new sample in the mc object

		// Fit
		localfit = opt.optimize (model, localdataset );
#ifdef DEBUG_BOOTSTRAP
		for (l=0; l<sample.size(); l++)
			std::cerr << " " << sample[l] << "\n";
		std::cerr << localfit[0] << " " << localfit[1] << " " << localfit[2] << "\n";
#endif

		// Get some characteristics of the localfit
		deviance = model->deviance ( localfit, localdataset );
		devianceresiduals = model->getDevianceResiduals ( localfit, localdataset );
		bootstrapsamples.setEst ( b, localfit, deviance );
		bootstrapsamples.setRpd ( b, model->getRpd( devianceresiduals, localfit, localdataset ) );
		bootstrapsamples.setRkd ( b, model->getRkd( devianceresiduals, localdataset ) );

		// Store what we need for the BCa stuff
		for (cut=0; cut<cuts.size(); cut++) {
			l_LF[cut][b] = model->leastfavourable ( localfit, localdataset, cuts[cut] );
#ifdef DEBUG_BOOTSTRAP
			if (l_LF[cut][b] != l_LF[cut][b]) {
				std::cerr << "deviance = " << deviance << "\n";
				std::cerr << "l_LF["<<cut<<"]["<<b<<"] = " << l_LF[cut][b] << "\n";
			}
#endif
			u_i[cut][b]  = model->getThres(localfit,cuts[cut]);
#ifdef DEBUG_BOOTSTRAP
			if (l_LF[cut][b]!= l_LF[cut][b]) {
				std::cerr << "u_i["<<cut<<"]["<<b<<"] = " << u_i[cut][b] << "\n";
			}
#endif
			bootstrapsamples.setThres(u_i[cut][b], b, cut);

			if (l_LF[cut][b] != l_LF[cut][b]) {
				// TODO: if l_LF is nan we don't take this sample
				// TODO: This is not the best solution but it works (kindof)
				b--;
				continue;
			}
		}
	}

	// Calculate BCa constants
	double bias, acc;
	for (cut=0; cut<cuts.size(); cut++) {
		determineBCa ( l_LF[cut], u_i[cut], initialthresholds[cut], &bias, &acc );
		bootstrapsamples.setBCa(cut, bias, acc );
	}

	delete localdataset;

	return bootstrapsamples;
}

JackKnifeList jackknifedata ( const PsiData * data, const PsiPsychometric* model )
{
	PsiOptimizer *opt = new PsiOptimizer( model, data );
	std::vector<double> mlestimate ( opt->optimize( model, data ) );
	delete opt;
	JackKnifeList jackknife ( data->getNblocks(), model->getNparams(), model->deviance(mlestimate, data) );
	PsiData * localdata;

	std::vector<double> x ( data->getNblocks()-1 );
	std::vector<int> k ( data->getNblocks()-1 );
	std::vector<int> n ( data->getNblocks()-1 );
	int i,j,l,exclude(0);

	for ( i=0; i<data->getNblocks(); i++ ) {
		exclude=i;
		j=0;
		for (l=0; l<data->getNblocks(); l++) {
			if (l!=exclude) {
				x[j] = data->getIntensity(l);
				k[j] = data->getNcorrect(l);
				n[j] = data->getNtrials(l);
				j++;
			}
		}

		localdata = new PsiData ( x,n,k,data->getNalternatives() );
		opt       = new PsiOptimizer ( model, localdata );

		mlestimate = opt->optimize( model, localdata );
		jackknife.setEst ( i, mlestimate, model->deviance(mlestimate,localdata) );

		delete localdata;
		delete opt;
	}

	return jackknife;
}

